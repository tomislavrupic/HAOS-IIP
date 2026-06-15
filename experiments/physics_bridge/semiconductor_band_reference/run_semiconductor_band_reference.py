#!/usr/bin/env python3
"""Semiconductor band-structure computational reference probe.

This sidecar runs a narrow toy band-structure reference check:

- a declared two-band direct-gap reference must reproduce its own band table;
- controls preserve some confounds while breaking gap, directness, curvature,
  dispersion, band identity, or k-ordering;
- real materials, carrier dynamics, transport, ab-initio band structure, and
  CST/HAOS derivation remain outside this probe.

It is not a physical semiconductor calculation and not a claim that HAOS-IIP
derives band structure.
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
CATALOG_PATH = ROOT / "band_catalog.json"
CSV_PATH = ROOT / "semiconductor_band_diagnostics.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "semiconductor_band_reference_report.md"


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
    "k_index",
    "k_reduced",
    "conduction_energy_ev",
    "valence_energy_ev",
    "reference_conduction_energy_ev",
    "reference_valence_energy_ev",
    "conduction_error_ev",
    "valence_error_ev",
    "band_gap_ev",
    "gap_relative_error",
    "conduction_edge_k",
    "valence_edge_k",
    "direct_gap_score",
    "curvature_score",
    "band_shape_score",
    "gap_score",
    "band_order_score",
    "symmetry_score",
    "semiconductor_reference_score",
    "status",
]


@dataclass(frozen=True)
class BandConfig:
    k_points: int = 81
    reference_gap_ev: float = 1.20
    conduction_curvature_scale_ev: float = 0.48
    valence_curvature_scale_ev: float = 0.36
    pass_max_reference_shape_error: float = 1.0e-12
    pass_max_gap_relative_error: float = 1.0e-12
    pass_min_direct_gap_score: float = 1.0
    pass_min_curvature_score: float = 0.999999
    pass_min_band_order_score: float = 1.0
    pass_min_symmetry_score: float = 0.999999
    pass_min_reference_score: float = 0.999999
    control_false_positive_score_gate: float = 0.95
    gap_score_error_scale: float = 0.05
    curvature_score_error_scale: float = 0.10
    band_shape_error_scale_ev: float = 0.08
    symmetry_error_scale_ev: float = 0.03
    edge_tolerance_k: float = 1.0e-12
    indirect_shift_k: float = 0.35
    scaled_gap_factor: float = 0.72
    wrong_curvature_factor: float = 0.35


@dataclass(frozen=True)
class BandPoint:
    k_index: int
    k_reduced: float
    reference_conduction_energy_ev: float
    reference_valence_energy_ev: float
    conduction_energy_ev: float
    valence_energy_ev: float


@dataclass(frozen=True)
class BandSummary:
    sample: str
    sample_role: str
    model_family: str
    point_count: int
    band_gap_ev: float
    gap_relative_error: float
    conduction_edge_k: float
    valence_edge_k: float
    direct_gap_score: float
    conduction_curvature: float
    valence_curvature_abs: float
    conduction_curvature_relative_error: float
    valence_curvature_relative_error: float
    curvature_score: float
    band_shape_rms_error_ev: float
    band_shape_score: float
    gap_score: float
    band_order_score: float
    symmetry_error_ev: float
    symmetry_score: float
    semiconductor_reference_score: float
    status: str


@dataclass(frozen=True)
class DiagnosticRow:
    sample: str
    sample_role: str
    model_family: str
    data_kind: str
    k_index: int | None
    k_reduced: float | None
    conduction_energy_ev: float | None
    valence_energy_ev: float | None
    reference_conduction_energy_ev: float | None
    reference_valence_energy_ev: float | None
    conduction_error_ev: float | None
    valence_error_ev: float | None
    band_gap_ev: float | None
    gap_relative_error: float | None
    conduction_edge_k: float | None
    valence_edge_k: float | None
    direct_gap_score: float | None
    curvature_score: float | None
    band_shape_score: float | None
    gap_score: float | None
    band_order_score: float | None
    symmetry_score: float | None
    semiconductor_reference_score: float | None
    status: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run semiconductor band-structure computational reference probe.")
    parser.add_argument("--k-points", type=int, default=BandConfig.k_points)
    parser.add_argument("--reference-gap-ev", type=float, default=BandConfig.reference_gap_ev)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated contract, catalog, CSV, status JSON, and markdown report.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> BandConfig:
    k_points = int(args.k_points)
    if k_points < 5 or k_points % 2 == 0:
        raise ValueError("k_points must be an odd integer >= 5 so k=0 is sampled.")
    return BandConfig(k_points=k_points, reference_gap_ev=float(args.reference_gap_ev))


def stable_json_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def rms(values: list[float]) -> float:
    if not values:
        return 0.0
    return math.sqrt(sum(value * value for value in values) / float(len(values)))


def k_grid(config: BandConfig) -> list[float]:
    half = config.k_points // 2
    return [float(index) / float(half) for index in range(-half, half + 1)]


def two_band_energy(k_value: float, gap_ev: float, conduction_scale: float, valence_scale: float, conduction_shift: float = 0.0, valence_shift: float = 0.0) -> tuple[float, float]:
    conduction = 0.5 * gap_ev + conduction_scale * (1.0 - math.cos(math.pi * (k_value - conduction_shift)))
    valence = -0.5 * gap_ev - valence_scale * (1.0 - math.cos(math.pi * (k_value - valence_shift)))
    return conduction, valence


def build_reference_catalog(config: BandConfig) -> list[BandPoint]:
    points: list[BandPoint] = []
    for index, k_value in enumerate(k_grid(config)):
        conduction, valence = two_band_energy(
            k_value,
            config.reference_gap_ev,
            config.conduction_curvature_scale_ev,
            config.valence_curvature_scale_ev,
        )
        points.append(
            BandPoint(
                k_index=index,
                k_reduced=k_value,
                reference_conduction_energy_ev=conduction,
                reference_valence_energy_ev=valence,
                conduction_energy_ev=conduction,
                valence_energy_ev=valence,
            )
        )
    return points


def replace_bands(reference: list[BandPoint], conduction_values: list[float], valence_values: list[float]) -> list[BandPoint]:
    return [
        BandPoint(
            k_index=point.k_index,
            k_reduced=point.k_reduced,
            reference_conduction_energy_ev=point.reference_conduction_energy_ev,
            reference_valence_energy_ev=point.reference_valence_energy_ev,
            conduction_energy_ev=float(conduction),
            valence_energy_ev=float(valence),
        )
        for point, conduction, valence in zip(reference, conduction_values, valence_values)
    ]


def build_scaled_gap_control(reference: list[BandPoint], config: BandConfig) -> list[BandPoint]:
    shift = 0.5 * (1.0 - config.scaled_gap_factor) * config.reference_gap_ev
    conduction = [point.reference_conduction_energy_ev - shift for point in reference]
    valence = [point.reference_valence_energy_ev + shift for point in reference]
    return replace_bands(reference, conduction, valence)


def build_metallic_overlap_control(reference: list[BandPoint], _config: BandConfig) -> list[BandPoint]:
    conduction = [point.reference_conduction_energy_ev - 0.68 for point in reference]
    valence = [point.reference_valence_energy_ev + 0.68 for point in reference]
    return replace_bands(reference, conduction, valence)


def build_indirect_gap_control(reference: list[BandPoint], config: BandConfig) -> list[BandPoint]:
    conduction: list[float] = []
    valence: list[float] = []
    for point in reference:
        c_value, v_value = two_band_energy(
            point.k_reduced,
            config.reference_gap_ev,
            config.conduction_curvature_scale_ev,
            config.valence_curvature_scale_ev,
            conduction_shift=config.indirect_shift_k,
        )
        conduction.append(c_value)
        valence.append(v_value)
    return replace_bands(reference, conduction, valence)


def build_wrong_curvature_control(reference: list[BandPoint], config: BandConfig) -> list[BandPoint]:
    conduction: list[float] = []
    valence: list[float] = []
    for point in reference:
        c_value, v_value = two_band_energy(
            point.k_reduced,
            config.reference_gap_ev,
            config.conduction_curvature_scale_ev * config.wrong_curvature_factor,
            config.valence_curvature_scale_ev * (2.0 - config.wrong_curvature_factor),
        )
        conduction.append(c_value)
        valence.append(v_value)
    return replace_bands(reference, conduction, valence)


def build_flat_band_control(reference: list[BandPoint], config: BandConfig) -> list[BandPoint]:
    conduction = [0.5 * config.reference_gap_ev for _point in reference]
    valence = [-0.5 * config.reference_gap_ev for _point in reference]
    return replace_bands(reference, conduction, valence)


def build_band_label_swap_control(reference: list[BandPoint], _config: BandConfig) -> list[BandPoint]:
    conduction = [point.reference_valence_energy_ev for point in reference]
    valence = [point.reference_conduction_energy_ev for point in reference]
    return replace_bands(reference, conduction, valence)


def build_k_order_shuffle_control(reference: list[BandPoint], _config: BandConfig) -> list[BandPoint]:
    conduction_values = [point.reference_conduction_energy_ev for point in reference]
    valence_values = [point.reference_valence_energy_ev for point in reference]
    even_values = conduction_values[::2]
    odd_values = conduction_values[1::2]
    shuffled_conduction = list(reversed(even_values)) + odd_values
    even_valence = valence_values[::2]
    odd_valence = valence_values[1::2]
    shuffled_valence = list(reversed(even_valence)) + odd_valence
    return replace_bands(reference, shuffled_conduction[: len(reference)], shuffled_valence[: len(reference)])


def edge_indices(points: list[BandPoint]) -> tuple[int, int]:
    conduction_index = min(range(len(points)), key=lambda idx: points[idx].conduction_energy_ev)
    valence_index = max(range(len(points)), key=lambda idx: points[idx].valence_energy_ev)
    return conduction_index, valence_index


def finite_curvature(values: list[float], edge_index: int, h_value: float) -> float:
    if edge_index <= 0 or edge_index >= len(values) - 1:
        return 0.0
    return (values[edge_index + 1] - 2.0 * values[edge_index] + values[edge_index - 1]) / max(h_value * h_value, 1.0e-300)


def symmetry_error(points: list[BandPoint]) -> float:
    by_k = {round(point.k_reduced, 12): point for point in points}
    errors: list[float] = []
    for point in points:
        opposite = by_k.get(round(-point.k_reduced, 12))
        if opposite is None:
            continue
        errors.append(abs(point.conduction_energy_ev - opposite.conduction_energy_ev))
        errors.append(abs(point.valence_energy_ev - opposite.valence_energy_ev))
    return sum(errors) / float(len(errors)) if errors else 0.0


def summarize_sample(
    sample: str,
    sample_role: str,
    model_family: str,
    points: list[BandPoint],
    reference: list[BandPoint],
    config: BandConfig,
) -> BandSummary:
    conduction_index, valence_index = edge_indices(points)
    ref_conduction_index, ref_valence_index = edge_indices(reference)
    conduction_edge_k = points[conduction_index].k_reduced
    valence_edge_k = points[valence_index].k_reduced
    reference_conduction_edge_k = reference[ref_conduction_index].k_reduced
    reference_valence_edge_k = reference[ref_valence_index].k_reduced
    band_gap = points[conduction_index].conduction_energy_ev - points[valence_index].valence_energy_ev
    gap_relative_error = abs(band_gap - config.reference_gap_ev) / max(abs(config.reference_gap_ev), 1.0e-300)
    direct_gap_score = (
        1.0
        if (
            abs(conduction_edge_k - valence_edge_k) <= config.edge_tolerance_k
            and abs(conduction_edge_k - reference_conduction_edge_k) <= config.edge_tolerance_k
            and abs(valence_edge_k - reference_valence_edge_k) <= config.edge_tolerance_k
        )
        else 0.0
    )

    h_value = abs(points[1].k_reduced - points[0].k_reduced)
    conduction_values = [point.conduction_energy_ev for point in points]
    valence_values = [point.valence_energy_ev for point in points]
    reference_conduction_values = [point.conduction_energy_ev for point in reference]
    reference_valence_values = [point.valence_energy_ev for point in reference]
    conduction_curvature = finite_curvature(conduction_values, conduction_index, h_value)
    valence_curvature_abs = abs(finite_curvature(valence_values, valence_index, h_value))
    ref_conduction_curvature = finite_curvature(reference_conduction_values, ref_conduction_index, h_value)
    ref_valence_curvature_abs = abs(finite_curvature(reference_valence_values, ref_valence_index, h_value))
    conduction_curvature_error = abs(conduction_curvature - ref_conduction_curvature) / max(abs(ref_conduction_curvature), 1.0e-300)
    valence_curvature_error = abs(valence_curvature_abs - ref_valence_curvature_abs) / max(abs(ref_valence_curvature_abs), 1.0e-300)
    curvature_error = rms([conduction_curvature_error, valence_curvature_error])

    shape_errors: list[float] = []
    for point in points:
        shape_errors.append(point.conduction_energy_ev - point.reference_conduction_energy_ev)
        shape_errors.append(point.valence_energy_ev - point.reference_valence_energy_ev)
    band_shape_rms = rms(shape_errors)

    gap_score = clamp01(1.0 - gap_relative_error / config.gap_score_error_scale)
    curvature_score = clamp01(1.0 - curvature_error / config.curvature_score_error_scale)
    band_shape_score = clamp01(1.0 - band_shape_rms / config.band_shape_error_scale_ev)
    ordered = sum(1 for point in points if point.conduction_energy_ev > point.valence_energy_ev)
    band_order_score = float(ordered) / float(len(points)) if points else 0.0
    sym_error = symmetry_error(points)
    symmetry_score = clamp01(1.0 - sym_error / config.symmetry_error_scale_ev)
    reference_score = min(
        gap_score,
        direct_gap_score,
        curvature_score,
        band_shape_score,
        band_order_score,
        symmetry_score,
    )

    if sample_role == "reference":
        status = (
            "REFERENCE_TOY_BAND_REPRODUCED"
            if (
                band_shape_rms <= config.pass_max_reference_shape_error
                and gap_relative_error <= config.pass_max_gap_relative_error
                and direct_gap_score >= config.pass_min_direct_gap_score
                and curvature_score >= config.pass_min_curvature_score
                and band_order_score >= config.pass_min_band_order_score
                and symmetry_score >= config.pass_min_symmetry_score
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

    return BandSummary(
        sample=sample,
        sample_role=sample_role,
        model_family=model_family,
        point_count=len(points),
        band_gap_ev=band_gap,
        gap_relative_error=gap_relative_error,
        conduction_edge_k=conduction_edge_k,
        valence_edge_k=valence_edge_k,
        direct_gap_score=direct_gap_score,
        conduction_curvature=conduction_curvature,
        valence_curvature_abs=valence_curvature_abs,
        conduction_curvature_relative_error=conduction_curvature_error,
        valence_curvature_relative_error=valence_curvature_error,
        curvature_score=curvature_score,
        band_shape_rms_error_ev=band_shape_rms,
        band_shape_score=band_shape_score,
        gap_score=gap_score,
        band_order_score=band_order_score,
        symmetry_error_ev=sym_error,
        symmetry_score=symmetry_score,
        semiconductor_reference_score=reference_score,
        status=status,
    )


def make_rows(
    sample: str,
    sample_role: str,
    model_family: str,
    points: list[BandPoint],
    summary: BandSummary,
) -> list[DiagnosticRow]:
    rows: list[DiagnosticRow] = [
        DiagnosticRow(
            sample=sample,
            sample_role=sample_role,
            model_family=model_family,
            data_kind="aggregate",
            k_index=None,
            k_reduced=None,
            conduction_energy_ev=None,
            valence_energy_ev=None,
            reference_conduction_energy_ev=None,
            reference_valence_energy_ev=None,
            conduction_error_ev=None,
            valence_error_ev=None,
            band_gap_ev=summary.band_gap_ev,
            gap_relative_error=summary.gap_relative_error,
            conduction_edge_k=summary.conduction_edge_k,
            valence_edge_k=summary.valence_edge_k,
            direct_gap_score=summary.direct_gap_score,
            curvature_score=summary.curvature_score,
            band_shape_score=summary.band_shape_score,
            gap_score=summary.gap_score,
            band_order_score=summary.band_order_score,
            symmetry_score=summary.symmetry_score,
            semiconductor_reference_score=summary.semiconductor_reference_score,
            status=summary.status,
        )
    ]
    for point in points:
        rows.append(
            DiagnosticRow(
                sample=sample,
                sample_role=sample_role,
                model_family=model_family,
                data_kind="band_point",
                k_index=point.k_index,
                k_reduced=point.k_reduced,
                conduction_energy_ev=point.conduction_energy_ev,
                valence_energy_ev=point.valence_energy_ev,
                reference_conduction_energy_ev=point.reference_conduction_energy_ev,
                reference_valence_energy_ev=point.reference_valence_energy_ev,
                conduction_error_ev=point.conduction_energy_ev - point.reference_conduction_energy_ev,
                valence_error_ev=point.valence_energy_ev - point.reference_valence_energy_ev,
                band_gap_ev=None,
                gap_relative_error=None,
                conduction_edge_k=None,
                valence_edge_k=None,
                direct_gap_score=None,
                curvature_score=None,
                band_shape_score=None,
                gap_score=None,
                band_order_score=None,
                symmetry_score=None,
                semiconductor_reference_score=None,
                status="BAND_POINT_RECORDED",
            )
        )
    return rows


def precommitment_contract(config: BandConfig) -> dict[str, Any]:
    return {
        "name": "Semiconductor band-structure computational reference probe",
        "status": "claim_gated_sidecar",
        "scope": {
            "implemented_fact": "Run a deterministic toy direct-gap two-band reference and component controls.",
            "design_choice": "Use declared cosine bands on a reduced k-grid with explicit gap and curvature parameters.",
            "heuristic": "Scores decompose gap, directness, curvature, band-shape, band-order, and k-symmetry checks.",
            "unverified_hypothesis": "No CST or HAOS semiconductor band-recovery hypothesis is tested here.",
        },
        "non_claims": [
            "not a physical semiconductor calculation",
            "not an ab-initio band-structure calculation",
            "not a carrier-transport or device model",
            "not a materials prediction",
            "not a CST derivation of band structure",
            "not evidence that HAOS-IIP recovers semiconductor physics",
        ],
        "reference_model": {
            "name": "declared toy direct-gap two-band reference",
            "k_domain": "reduced unitless grid [-1, 1]",
            "formula": {
                "conduction": "E_c(k) = Eg/2 + A_c * (1 - cos(pi*k))",
                "valence": "E_v(k) = -Eg/2 - A_v * (1 - cos(pi*k))",
            },
            "parameters": {
                "k_points": config.k_points,
                "reference_gap_ev": config.reference_gap_ev,
                "conduction_curvature_scale_ev": config.conduction_curvature_scale_ev,
                "valence_curvature_scale_ev": config.valence_curvature_scale_ev,
            },
            "note": "This declared model is the calibration target; it is not derived from HAOS-IIP.",
        },
        "controls": {
            "scaled_gap_control": {
                "preserve": ["curvature shape", "direct gap location", "k symmetry", "band ordering"],
                "destroy": ["absolute band gap"],
                "expected_response": "gap_score and band_shape_score decrease while curvature_score remains high",
            },
            "metallic_overlap_control": {
                "preserve": ["curvature shape", "k symmetry"],
                "destroy": ["positive band gap", "band ordering near the edge"],
                "expected_response": "gap_score and band_order_score decrease",
            },
            "indirect_gap_control": {
                "preserve": ["positive gap", "curvature scales", "band ordering"],
                "destroy": ["direct gap alignment", "k symmetry"],
                "expected_response": "direct_gap_score and symmetry_score decrease",
            },
            "wrong_curvature_control": {
                "preserve": ["gap", "direct gap location", "band ordering", "k symmetry"],
                "destroy": ["effective-mass/curvature proxy"],
                "expected_response": "curvature_score and band_shape_score decrease",
            },
            "flat_band_control": {
                "preserve": ["gap", "direct gap location", "band ordering", "k symmetry"],
                "destroy": ["dispersion and curvature"],
                "expected_response": "curvature_score and band_shape_score decrease",
            },
            "band_label_swap_control": {
                "preserve": ["energy magnitudes", "k symmetry"],
                "destroy": ["conduction/valence identity", "band ordering"],
                "expected_response": "band_order_score and gap_score decrease",
            },
            "k_order_shuffle_control": {
                "preserve": ["energy value marginals"],
                "destroy": ["smooth k-ordering", "edge identity", "curvature"],
                "expected_response": "band_shape_score, curvature_score, and symmetry_score decrease",
            },
        },
        "control_parameters": {
            "indirect_shift_k": config.indirect_shift_k,
            "scaled_gap_factor": config.scaled_gap_factor,
            "wrong_curvature_factor": config.wrong_curvature_factor,
        },
        "thresholds": {
            "pass_max_reference_shape_error": config.pass_max_reference_shape_error,
            "pass_max_gap_relative_error": config.pass_max_gap_relative_error,
            "pass_min_direct_gap_score": config.pass_min_direct_gap_score,
            "pass_min_curvature_score": config.pass_min_curvature_score,
            "pass_min_band_order_score": config.pass_min_band_order_score,
            "pass_min_symmetry_score": config.pass_min_symmetry_score,
            "pass_min_reference_score": config.pass_min_reference_score,
            "control_false_positive_score_gate": config.control_false_positive_score_gate,
            "gap_score_error_scale": config.gap_score_error_scale,
            "curvature_score_error_scale": config.curvature_score_error_scale,
            "band_shape_error_scale_ev": config.band_shape_error_scale_ev,
            "symmetry_error_scale_ev": config.symmetry_error_scale_ev,
        },
        "verdict_logic": {
            "reference_pass": [
                "reference band-shape RMS error <= pass_max_reference_shape_error",
                "gap relative error <= pass_max_gap_relative_error",
                "direct_gap_score >= pass_min_direct_gap_score",
                "curvature_score >= pass_min_curvature_score",
                "band_order_score >= pass_min_band_order_score",
                "symmetry_score >= pass_min_symmetry_score",
                "all controls score < control_false_positive_score_gate",
            ],
            "always_report": [
                "CST_SEMICONDUCTOR_STATUS_OPEN",
                "HAOS_DERIVATION_NOT_TESTED",
                "TOY_BAND_STRUCTURE_ONLY",
            ],
        },
        "future_candidate_rejection_conditions": [
            "candidate imports the declared cosine band formula by construction",
            "candidate imports Eg or curvature parameters by construction",
            "candidate passes only by fitting after seeing the reference table",
            "candidate uses band labels or k ordering unavailable to controls",
            "candidate claims real material, transport, effective-mass, or device recovery from this toy reference",
        ],
    }


def sample_definitions(reference: list[BandPoint], config: BandConfig) -> list[tuple[str, str, str, list[BandPoint]]]:
    return [
        ("toy_direct_gap_reference", "reference", "declared_two_band_reference", reference),
        ("scaled_gap_control", "negative_control", "gap_scale_control", build_scaled_gap_control(reference, config)),
        ("metallic_overlap_control", "negative_control", "gap_destroying_control", build_metallic_overlap_control(reference, config)),
        ("indirect_gap_control", "negative_control", "directness_destroying_control", build_indirect_gap_control(reference, config)),
        ("wrong_curvature_control", "negative_control", "curvature_destroying_control", build_wrong_curvature_control(reference, config)),
        ("flat_band_control", "negative_control", "dispersion_destroying_control", build_flat_band_control(reference, config)),
        ("band_label_swap_control", "negative_control", "identity_destroying_control", build_band_label_swap_control(reference, config)),
        ("k_order_shuffle_control", "negative_control", "k_order_destroying_control", build_k_order_shuffle_control(reference, config)),
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


def run_probe(config: BandConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = precommitment_contract(config)
    contract["contract_hash"] = stable_json_hash("semiconductor_contract", contract)
    write_json(output_dir / CONTRACT_PATH.name, contract)

    reference = build_reference_catalog(config)
    write_json(output_dir / CATALOG_PATH.name, [asdict(point) for point in reference])

    summaries: list[BandSummary] = []
    rows: list[DiagnosticRow] = []
    for sample, sample_role, model_family, points in sample_definitions(reference, config):
        summary = summarize_sample(sample, sample_role, model_family, points, reference, config)
        summaries.append(summary)
        rows.extend(make_rows(sample, sample_role, model_family, points, summary))
    write_csv(output_dir / CSV_PATH.name, rows)

    reference_summary = next(summary for summary in summaries if summary.sample_role == "reference")
    control_false_positive_count = sum(
        1 for summary in summaries if summary.sample_role == "negative_control" and summary.status == "CONTROL_FALSE_POSITIVE"
    )
    reference_pass = (
        reference_summary.status == "REFERENCE_TOY_BAND_REPRODUCED"
        and control_false_positive_count == 0
    )
    labels = [
        "SEMICONDUCTOR_REFERENCE_SANITY_PASS" if reference_pass else "SEMICONDUCTOR_REFERENCE_SANITY_FAIL",
        "CST_SEMICONDUCTOR_STATUS_OPEN",
        "HAOS_DERIVATION_NOT_TESTED",
        "TOY_BAND_STRUCTURE_ONLY",
        "NO_PHYSICAL_EXPERIMENT",
    ]
    result: dict[str, Any] = {
        "labels": labels,
        "contract_hash": contract["contract_hash"],
        "config": asdict(config),
        "reference_model": contract["reference_model"],
        "summaries": [asdict(summary) for summary in summaries],
        "control_false_positive_count": control_false_positive_count,
        "cst_mechanism_slot": "EMPTY_NOT_TESTED",
        "real_material_status": "OPEN_NOT_MODELED",
        "outputs": {
            "precommitment_contract": repo_rel(output_dir / CONTRACT_PATH.name),
            "band_catalog": repo_rel(output_dir / CATALOG_PATH.name),
            "semiconductor_band_diagnostics": repo_rel(output_dir / CSV_PATH.name),
            "bridge_status": repo_rel(output_dir / STATUS_PATH.name),
            "semiconductor_band_reference_report": repo_rel(output_dir / REPORT_PATH.name),
        },
        "limitations": contract["non_claims"],
    }
    result["result_hash"] = stable_json_hash("semiconductor_band_result", result)
    write_json(output_dir / STATUS_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result, summaries)
    return result


def write_report(path: Path, result: dict[str, Any], summaries: list[BandSummary]) -> None:
    reference = next(summary for summary in summaries if summary.sample_role == "reference")
    controls = [summary for summary in summaries if summary.sample_role == "negative_control"]
    lines = [
        "# Semiconductor Band-Structure Computational Reference Probe",
        "",
        "Implemented fact: this sidecar runs a deterministic toy direct-gap two-band reference and component controls.",
        "Design choice: the reference uses declared cosine bands on a reduced unitless k-grid.",
        "Heuristic: the reference score is the minimum of gap, directness, curvature, band-shape, band-order, and symmetry scores.",
        "Unverified hypothesis: no CST or HAOS semiconductor band-recovery hypothesis is tested here.",
        "",
        "## Verdict Labels",
    ]
    lines.extend(f"- {label}" for label in result["labels"])
    lines.extend(
        [
            "",
            "## Reference",
            f"- reference gap: `{result['config']['reference_gap_ev']:.12f} eV`",
            f"- k-point count: `{reference.point_count}`",
            f"- measured band gap: `{reference.band_gap_ev:.12f} eV`",
            f"- conduction edge k: `{reference.conduction_edge_k:.12f}`",
            f"- valence edge k: `{reference.valence_edge_k:.12f}`",
            f"- reference score: `{reference.semiconductor_reference_score:.12f}`",
            "",
            "## Controls",
        ]
    )
    for control in controls:
        lines.append(
            "- {sample}: score `{score:.12f}`, gap `{gap:.12f}`, direct `{direct:.12f}`, curvature `{curvature:.12f}`, shape `{shape:.12f}`, order `{order:.12f}`, symmetry `{symmetry:.12f}`, status `{status}`".format(
                sample=control.sample,
                score=control.semiconductor_reference_score,
                gap=control.gap_score,
                direct=control.direct_gap_score,
                curvature=control.curvature_score,
                shape=control.band_shape_score,
                order=control.band_order_score,
                symmetry=control.symmetry_score,
                status=control.status,
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "- This is not a physical semiconductor calculation.",
            "- This does not model real materials, doping, defects, many-body effects, transport, devices, or carrier dynamics.",
            "- This does not derive band structure from CST or HAOS-IIP.",
            "- This does not change the frozen CST v0.2.2 checkpoint.",
            "- `SEMICONDUCTOR_REFERENCE_SANITY_PASS` means only that the toy reference harness behaves as expected.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_probe(make_config(args), args.output_dir)
    print(json.dumps({"labels": result["labels"], "result_hash": result["result_hash"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
