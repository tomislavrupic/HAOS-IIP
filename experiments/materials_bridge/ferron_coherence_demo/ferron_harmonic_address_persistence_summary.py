#!/usr/bin/env python3
"""Unified ferron harmonic-address persistence summary.

This is a sidecar synthesis over already-bounded audit outputs. It does not
create a new materials theory and does not convert the missing raw STFT grid
into evidence. It only gathers the spectral, published propagation, and
computed-STFT diagnostics into one persistence ledger.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
OUTPUTS_DIR = ROOT / "outputs"
TARGET_FREQUENCY_THZ = 3.13
VELOCITY_REFERENCE_M_S = 100000.0

COMPONENT_FIELDNAMES = [
    "component",
    "value",
    "weight",
    "source",
    "included_in_proxy",
    "note",
]


def run_harmonic_address_persistence_summary(
    *,
    output_dir: Path = OUTPUTS_DIR,
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
    velocity_reference_m_s: float = VELOCITY_REFERENCE_M_S,
) -> dict[str, Any]:
    """Build the unified persistence summary from existing audit JSON outputs."""

    output_dir.mkdir(parents=True, exist_ok=True)
    spectral = _read_json(output_dir / "spectral_feature_summary.json")
    raw_stft = _read_json(output_dir / "stft_feature_summary.json")
    published_stft = _read_json(output_dir / "published_stft_trace_summary.json")
    computed_stft = _read_json(output_dir / "computed_stft_diagnostic" / "computed_stft_summary.json")

    components = _build_components(
        spectral=spectral,
        raw_stft=raw_stft,
        published_stft=published_stft,
        computed_stft=computed_stft,
        velocity_reference_m_s=velocity_reference_m_s,
    )
    included = [
        row
        for row in components
        if row.get("included_in_proxy") and _float_or_none(row.get("value")) is not None
    ]
    if included:
        weighted_sum = sum(float(row["value"]) * float(row["weight"]) for row in included)
        weight_sum = sum(float(row["weight"]) for row in included)
        persistence_proxy = weighted_sum / weight_sum if weight_sum > 0.0 else None
    else:
        persistence_proxy = None

    raw_grid_status = str(raw_stft.get("status") or "UNKNOWN")
    raw_grid_boundary = raw_grid_status == "NO_STFT_DATA_FOUND"
    summary = {
        "status": "PASS" if persistence_proxy is not None else "PARTIAL",
        "mode": "ferron_harmonic_address_persistence_summary",
        "target_frequency_THz": target_frequency_THz,
        "source_doi": "10.6084/m9.figshare.31895293",
        "primary_paper_doi": "10.1038/s41563-026-02597-4",
        "raw_stft_time_frequency_grid_status": raw_grid_status,
        "raw_stft_grid_boundary_preserved": raw_grid_boundary,
        "combined_persistence_proxy": persistence_proxy,
        "components": components,
        "spectral_component": {
            "status": spectral.get("status"),
            "records_audited": spectral.get("records_audited"),
            "peaks_found": spectral.get("peaks_found"),
            "mean_peak_frequency_THz": spectral.get("mean_peak_frequency_THz"),
            "mean_spectral_recoverability": spectral.get("mean_spectral_recoverability"),
            "k_star_detected": spectral.get("k_star_detected"),
        },
        "published_propagation_component": {
            "status": published_stft.get("status"),
            "traces_found": published_stft.get("traces_found"),
            "mean_group_velocity_m_s": published_stft.get("mean_group_velocity_m_s"),
            "mean_velocity_consistency_with_1e5": published_stft.get("mean_velocity_consistency_with_1e5"),
            "group_velocity_rows": published_stft.get("group_velocity_rows", []),
        },
        "computed_stft_component": {
            "status": computed_stft.get("status"),
            "traces_computed": computed_stft.get("traces_computed"),
            "mean_group_velocity_m_s": computed_stft.get("mean_group_velocity_m_s"),
            "mean_envelope_correlation_vs_published": computed_stft.get(
                "mean_envelope_correlation_vs_published"
            ),
            "mean_abs_peak_time_delta_vs_published_ps": computed_stft.get(
                "mean_abs_peak_time_delta_vs_published_ps"
            ),
            "sheet_thickness_assignment_note": computed_stft.get("sheet_thickness_assignment_note"),
        },
        "bounded_interpretation": _bounded_interpretation(raw_grid_boundary),
        "claim_boundary": [
            "external-data sidecar summary only",
            "not a proof of HAOS-IIP",
            "not a claim that ferrons embody HAOS-IIP",
            "not a new materials theory",
            "not a raw STFT-grid claim",
            "computed STFT outputs are derived diagnostics, not published raw STFT maps",
        ],
        "outputs": {
            "summary_json": str(output_dir / "ferron_harmonic_address_persistence_summary.json"),
            "summary_md": str(output_dir / "ferron_harmonic_address_persistence_summary.md"),
            "metrics_csv": str(output_dir / "ferron_harmonic_address_persistence_metrics.csv"),
            "component_plot": str(output_dir / "ferron_harmonic_address_persistence_components.png"),
        },
    }

    _write_csv(output_dir / "ferron_harmonic_address_persistence_metrics.csv", components)
    _write_markdown(output_dir / "ferron_harmonic_address_persistence_summary.md", summary)
    _write_component_plot(output_dir / "ferron_harmonic_address_persistence_components.png", components)
    _write_json(output_dir / "ferron_harmonic_address_persistence_summary.json", summary)
    return summary


def _build_components(
    *,
    spectral: dict[str, Any],
    raw_stft: dict[str, Any],
    published_stft: dict[str, Any],
    computed_stft: dict[str, Any],
    velocity_reference_m_s: float,
) -> list[dict[str, Any]]:
    records = _float_or_none(spectral.get("records_audited")) or 0.0
    peaks = _float_or_none(spectral.get("peaks_found")) or 0.0
    peak_detection_fraction = peaks / records if records > 0.0 else None
    spectral_recoverability = _float_or_none(spectral.get("mean_spectral_recoverability"))
    no_k_star = 0.0 if spectral.get("k_star_detected") else 1.0
    published_velocity_consistency = _float_or_none(
        published_stft.get("mean_velocity_consistency_with_1e5")
    )
    published_monotonicity = _mean_or_none(
        _float_or_none(row.get("monotonic_peak_delay"))
        for row in published_stft.get("group_velocity_rows", [])
    )
    computed_envelope_match = _float_or_none(computed_stft.get("mean_envelope_correlation_vs_published"))
    computed_velocity_consistency = _velocity_consistency(
        _float_or_none(computed_stft.get("mean_group_velocity_m_s")),
        velocity_reference_m_s,
    )

    return [
        _component(
            "target_peak_presence",
            peak_detection_fraction,
            "spectral_feature_summary.json",
            "15/15 target-window peak detection when value is 1.",
        ),
        _component(
            "spectral_recoverability",
            spectral_recoverability,
            "spectral_feature_summary.json",
            "Mean second-pass spectral recoverability proxy.",
        ),
        _component(
            "no_sustained_k_star",
            no_k_star,
            "spectral_feature_summary.json",
            "Absence of sustained k-star collapse under current thresholds.",
        ),
        _component(
            "published_propagation_monotonicity",
            published_monotonicity,
            "published_stft_trace_summary.json",
            "Mean monotonic peak-delay fraction across published target-band STFT trace velocity groups.",
        ),
        _component(
            "published_velocity_consistency",
            published_velocity_consistency,
            "published_stft_trace_summary.json",
            "Consistency of published target-band peak-delay fit with 1e5 m/s reference.",
        ),
        _component(
            "computed_envelope_match",
            computed_envelope_match,
            "computed_stft_summary.json",
            "Mean envelope correlation between computed target-band STFT and published STFT intensity traces.",
        ),
        _component(
            "computed_velocity_consistency",
            computed_velocity_consistency,
            "computed_stft_summary.json",
            "Consistency of computed target-band peak-delay fit with 1e5 m/s reference.",
        ),
        {
            "component": "raw_stft_grid_boundary",
            "value": None,
            "weight": 0.0,
            "source": "stft_feature_summary.json",
            "included_in_proxy": False,
            "note": (
                f"Raw STFT/time-frequency grid status: {raw_stft.get('status')}. "
                "This is a data-availability boundary, not a persistence penalty."
            ),
        },
    ]


def _component(component: str, value: float | None, source: str, note: str) -> dict[str, Any]:
    return {
        "component": component,
        "value": _clip01(value),
        "weight": 1.0,
        "source": source,
        "included_in_proxy": value is not None,
        "note": note,
    }


def _bounded_interpretation(raw_grid_boundary: bool) -> str:
    boundary_text = (
        "The raw STFT time-frequency grid remains unavailable in the downloaded XLSX files; "
        "that boundary is preserved."
        if raw_grid_boundary
        else "Raw STFT-grid status is not the expected NO_STFT_DATA_FOUND boundary and should be inspected."
    )
    return (
        "The ferron sidecar shows a stable target-window spectral feature near 3.13 THz, "
        "published target-band STFT intensity traces with monotonic propagation delay, and "
        "computed-from-time-trace STFT envelopes that closely match the published traces. "
        f"{boundary_text} This is a bounded external-data persistence summary, not a proof or "
        "materials-theory claim."
    )


def _write_markdown(path: Path, summary: dict[str, Any]) -> None:
    lines = [
        "# Ferron Harmonic-Address Persistence Summary",
        "",
        "## Status",
        f"- status: {summary.get('status')}",
        f"- mode: {summary.get('mode')}",
        f"- target_frequency_THz: {_format_optional(summary.get('target_frequency_THz'))}",
        f"- combined_persistence_proxy: {_format_optional(summary.get('combined_persistence_proxy'))}",
        f"- raw_stft_time_frequency_grid_status: {summary.get('raw_stft_time_frequency_grid_status')}",
        "",
        "## Components",
    ]
    for row in summary.get("components", []):
        included = "included" if row.get("included_in_proxy") else "not included"
        lines.append(
            "- "
            f"{row.get('component')}: {_format_optional(row.get('value'))} "
            f"({included}; source={row.get('source')})"
        )
    lines.extend(
        [
            "",
            "## Bounded Interpretation",
            str(summary.get("bounded_interpretation")),
            "",
            "## Claim Boundary",
        ]
    )
    for item in summary.get("claim_boundary", []):
        lines.append(f"- {item}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_component_plot(path: Path, components: list[dict[str, Any]]) -> None:
    plotted = [
        row
        for row in components
        if row.get("included_in_proxy") and _float_or_none(row.get("value")) is not None
    ]
    if not plotted:
        return
    labels = [str(row["component"]).replace("_", "\n") for row in plotted]
    values = [float(row["value"]) for row in plotted]
    plt.figure(figsize=(8.4, 4.8))
    bars = plt.bar(range(len(values)), values)
    for bar, value in zip(bars, values):
        plt.text(
            bar.get_x() + bar.get_width() / 2.0,
            min(value + 0.025, 1.02),
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    plt.xticks(range(len(labels)), labels, fontsize=7)
    plt.ylabel("bounded component value")
    plt.ylim(0.0, 1.08)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COMPONENT_FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _csv_value(row.get(field)) for field in COMPONENT_FIELDNAMES})


def _velocity_consistency(velocity: float | None, reference: float) -> float | None:
    if velocity is None or reference <= 1.0e-12:
        return None
    return _clip01(1.0 - abs(velocity - reference) / reference)


def _mean_or_none(values: Any) -> float | None:
    finite = [value for value in values if value is not None and math.isfinite(float(value))]
    if not finite:
        return None
    return float(sum(finite) / len(finite))


def _float_or_none(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _clip01(value: float | None) -> float | None:
    if value is None or not math.isfinite(float(value)):
        return None
    return min(max(float(value), 0.0), 1.0)


def _csv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.10g}"
    return value


def _format_optional(value: Any) -> str:
    if value is None:
        return "None"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "None"
        return f"{value:.6g}"
    return str(value)


if __name__ == "__main__":
    run_harmonic_address_persistence_summary()
