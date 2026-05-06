#!/usr/bin/env python3
"""Second-pass spectral feature audit for the ferron materials bridge.

This module inspects parsed ferron / NbOI2 spectra for the target narrow-band
feature near 3.13 THz. It reports peak stability, amplitude retention,
linewidth proxies, and peak-to-background clarity. It does not modify the
first-pass k-star result and does not force collapse.
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
import numpy as np


TARGET_FREQUENCY_THZ = 3.13
SEARCH_WINDOW_THZ = 0.20
COLLAPSE_THRESHOLD = 0.70
SUSTAIN_STEPS = 2

SUMMARY_FIELDNAMES = [
    "source_file",
    "sheet_name",
    "sample_id",
    "condition_label",
    "distance_um",
    "temperature_K",
    "thickness_nm",
    "fluence_mJ_cm2",
    "target_frequency_THz",
    "peak_found",
    "peak_frequency_THz",
    "peak_amplitude",
    "linewidth_proxy_THz",
    "peak_to_background",
    "frequency_drift_from_target_THz",
    "amplitude_retention",
    "frequency_stability",
    "linewidth_stability",
    "spectral_recoverability",
]

GROUP_FIELDNAMES = [
    "group_key",
    "n_records",
    "mean_peak_frequency_THz",
    "std_peak_frequency_THz",
    "mean_peak_amplitude",
    "mean_linewidth_proxy_THz",
    "mean_peak_to_background",
    "mean_spectral_recoverability",
    "min_spectral_recoverability",
    "max_spectral_recoverability",
]


def detect_target_peak(
    frequency: np.ndarray,
    amplitude: np.ndarray,
    *,
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
    search_window_THz: float = SEARCH_WINDOW_THZ,
) -> dict[str, Any]:
    """Detect the strongest amplitude inside the target frequency window."""

    freq = np.asarray(frequency, dtype=float)
    values = np.abs(np.asarray(amplitude, dtype=float))
    finite = np.isfinite(freq) & np.isfinite(values)
    freq = freq[finite]
    values = values[finite]
    if freq.size == 0 or values.size == 0 or freq.size != values.size:
        return _missing_peak("empty_or_invalid_spectrum")

    window = np.abs(freq - target_frequency_THz) <= search_window_THz
    if not np.any(window):
        return _missing_peak("no_frequency_points_in_target_window")

    window_indices = np.where(window)[0]
    local_peak = int(np.nanargmax(values[window]))
    peak_index = int(window_indices[local_peak])
    return {
        "peak_frequency_THz": float(freq[peak_index]),
        "peak_amplitude": float(values[peak_index]),
        "peak_index": peak_index,
        "peak_found": True,
        "reason": None,
    }


def estimate_linewidth_proxy(
    frequency: np.ndarray,
    amplitude: np.ndarray,
    peak_index: int | None,
    *,
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
    search_window_THz: float = SEARCH_WINDOW_THZ,
) -> tuple[float | None, str | None]:
    """Estimate a deterministic FWHM-like width around the detected peak."""

    if peak_index is None:
        return None, "peak_index_missing"

    freq = np.asarray(frequency, dtype=float)
    values = np.abs(np.asarray(amplitude, dtype=float))
    finite = np.isfinite(freq) & np.isfinite(values)
    original_indices = np.where(finite)[0]
    freq = freq[finite]
    values = values[finite]
    if freq.size < 4 or values.size != freq.size:
        return None, "spectrum_too_short_for_linewidth"

    mapped = np.where(original_indices == peak_index)[0]
    if mapped.size == 0:
        return None, "peak_index_not_in_finite_spectrum"
    local_peak_index = int(mapped[0])

    window = np.abs(freq - target_frequency_THz) <= search_window_THz
    if not np.any(window):
        return None, "no_target_window_for_linewidth"

    outside = ~window
    if np.any(outside):
        background = float(np.nanmedian(values[outside]))
    else:
        background = float(np.nanmin(values[window]))
    adjusted = np.maximum(values - background, 0.0)
    peak_height = float(adjusted[local_peak_index])
    if not math.isfinite(peak_height) or peak_height <= 1.0e-12:
        return None, "peak_height_not_above_background"

    half_height = 0.5 * peak_height
    left_cross = _half_max_crossing(freq, adjusted, local_peak_index, half_height, direction=-1)
    right_cross = _half_max_crossing(freq, adjusted, local_peak_index, half_height, direction=1)
    if left_cross is None:
        return None, "left_half_max_crossing_not_found"
    if right_cross is None:
        return None, "right_half_max_crossing_not_found"
    width = float(right_cross - left_cross)
    if not math.isfinite(width) or width <= 0.0:
        return None, "invalid_linewidth_width"
    return width, None


def estimate_peak_to_background(
    frequency: np.ndarray,
    amplitude: np.ndarray,
    peak_amplitude: float | None,
    *,
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
    search_window_THz: float = SEARCH_WINDOW_THZ,
) -> tuple[float | None, str | None]:
    """Return peak amplitude divided by median amplitude outside the target window."""

    if peak_amplitude is None or not math.isfinite(float(peak_amplitude)):
        return None, "peak_amplitude_missing"
    freq = np.asarray(frequency, dtype=float)
    values = np.abs(np.asarray(amplitude, dtype=float))
    finite = np.isfinite(freq) & np.isfinite(values)
    freq = freq[finite]
    values = values[finite]
    if freq.size == 0:
        return None, "empty_spectrum_for_background"
    outside = np.abs(freq - target_frequency_THz) > search_window_THz
    if not np.any(outside):
        return None, "no_background_points_outside_target_window"
    background = float(np.nanmedian(values[outside]))
    if not math.isfinite(background) or background <= 1.0e-12:
        return None, "background_median_zero_or_invalid"
    return float(peak_amplitude) / background, None


def compute_spectral_recoverability(
    *,
    amplitude_retention: float | None,
    frequency_stability: float | None,
    linewidth_stability: float | None,
    peak_to_background_stability: float | None,
) -> float | None:
    """Average available bounded components into a spectral recoverability proxy.

    Formula:

    spectral_recoverability = mean([
        clipped amplitude retention,
        clipped peak-frequency stability,
        clipped linewidth stability when available,
        clipped peak-to-background stability when available,
    ])

    Missing linewidth or background components are not fabricated; they are
    omitted and recorded in the JSON notes.
    """

    components = [
        _clip01(amplitude_retention),
        _clip01(frequency_stability),
        _clip01(linewidth_stability),
        _clip01(peak_to_background_stability),
    ]
    usable = [value for value in components if value is not None]
    if not usable:
        return None
    return float(sum(usable) / len(usable))


def audit_distance_decay(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Summarize amplitude retention versus distance where distance metadata exists."""

    by_group: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if _float_or_none(row.get("distance_um")) is None:
            continue
        by_group.setdefault(str(row.get("_baseline_group_key", "unknown")), []).append(row)

    summaries: list[dict[str, Any]] = []
    for group_key, group_rows in sorted(by_group.items()):
        ordered = sorted(group_rows, key=lambda row: (_float_or_none(row.get("distance_um")) or 0.0))
        if len(ordered) < 2:
            continue
        amplitudes = [_float_or_none(row.get("peak_amplitude")) for row in ordered]
        distances = [_float_or_none(row.get("distance_um")) for row in ordered]
        valid = [(d, a) for d, a in zip(distances, amplitudes) if d is not None and a is not None]
        if len(valid) < 2:
            continue
        nonincreasing = 0
        for (_, left_amp), (_, right_amp) in zip(valid, valid[1:]):
            if right_amp <= left_amp + 1.0e-12:
                nonincreasing += 1
        consistency = nonincreasing / max(len(valid) - 1, 1)
        baseline = max(valid[0][1], 1.0e-12)
        summaries.append(
            {
                "group_key": group_key,
                "distance_min_um": valid[0][0],
                "distance_max_um": valid[-1][0],
                "amplitude_retention_at_max_distance": min(max(valid[-1][1] / baseline, 0.0), 1.0),
                "raw_amplitude_ratio_at_max_distance": valid[-1][1] / baseline,
                "monotonic_decay_consistency": consistency,
                "n_distance_records": len(valid),
            }
        )
    return summaries


def audit_condition_groups(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group audited records by available source, sheet, sample, and condition metadata."""

    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        key = _condition_group_key(row)
        grouped.setdefault(key, []).append(row)

    group_rows: list[dict[str, Any]] = []
    for key, records in sorted(grouped.items()):
        peak_freq = _finite_values(record.get("peak_frequency_THz") for record in records)
        peak_amp = _finite_values(record.get("peak_amplitude") for record in records)
        linewidth = _finite_values(record.get("linewidth_proxy_THz") for record in records)
        peak_to_background = _finite_values(record.get("peak_to_background") for record in records)
        recoverability = _finite_values(record.get("spectral_recoverability") for record in records)
        group_rows.append(
            {
                "group_key": key,
                "n_records": len(records),
                "mean_peak_frequency_THz": _mean_or_none(peak_freq),
                "std_peak_frequency_THz": _std_or_none(peak_freq),
                "mean_peak_amplitude": _mean_or_none(peak_amp),
                "mean_linewidth_proxy_THz": _mean_or_none(linewidth),
                "mean_peak_to_background": _mean_or_none(peak_to_background),
                "mean_spectral_recoverability": _mean_or_none(recoverability),
                "min_spectral_recoverability": min(recoverability) if recoverability else None,
                "max_spectral_recoverability": max(recoverability) if recoverability else None,
            }
        )
    return group_rows


def summarize_standout_features(
    rows: list[dict[str, Any]],
    *,
    k_star_level: Any | None,
    distance_decay: list[dict[str, Any]],
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
) -> tuple[list[str], str]:
    found = [row for row in rows if row.get("peak_found")]
    if not found:
        return [], "No usable spectral target peak was detected in the parsed records."

    peak_frequencies = _finite_values(row.get("peak_frequency_THz") for row in found)
    recoverability = _finite_values(row.get("spectral_recoverability") for row in found)
    linewidth = _finite_values(row.get("linewidth_proxy_THz") for row in found)
    peak_to_background = _finite_values(row.get("peak_to_background") for row in found)
    drift = _finite_values(row.get("frequency_drift_from_target_THz") for row in found)

    features = [
        (
            "Target-window peak detected in "
            f"{len(found)}/{len(rows)} audited spectral records near {target_frequency_THz:.2f} THz."
        )
    ]
    if peak_frequencies:
        features.append(
            "Detected peak-frequency range is "
            f"{min(peak_frequencies):.6g}-{max(peak_frequencies):.6g} THz "
            f"(mean {sum(peak_frequencies) / len(peak_frequencies):.6g} THz)."
        )
    if drift:
        features.append(
            "Absolute drift from the 3.13 THz target ranges from "
            f"{min(abs(value) for value in drift):.6g} to {max(abs(value) for value in drift):.6g} THz."
        )
    if linewidth:
        features.append(
            "FWHM-like linewidth proxy is available for "
            f"{len(linewidth)} records, with mean {sum(linewidth) / len(linewidth):.6g} THz."
        )
    else:
        features.append("Linewidth proxy could not be estimated for the audited spectra.")
    if peak_to_background:
        features.append(
            "Peak-to-background ratio is available for "
            f"{len(peak_to_background)} records, with mean {sum(peak_to_background) / len(peak_to_background):.6g}."
        )
    if distance_decay:
        consistency_values = _finite_values(item.get("monotonic_decay_consistency") for item in distance_decay)
        retention_values = _finite_values(item.get("amplitude_retention_at_max_distance") for item in distance_decay)
        if consistency_values and retention_values:
            features.append(
                "Distance-series target peaks remain detectable; amplitude retention at maximum parsed distance is "
                f"{min(retention_values):.6g}-{max(retention_values):.6g}, while monotonic decay consistency is "
                f"{min(consistency_values):.6g}-{max(consistency_values):.6g}."
            )
    if recoverability:
        features.append(
            "Spectral recoverability ranges from "
            f"{min(recoverability):.6g} to {max(recoverability):.6g} across audited records."
        )

    if k_star_level is None:
        interpretation = (
            "Standout spectral feature: a target-window peak near 3.13 THz is detected in the parsed "
            "spectral records. Peak frequency remains stable under the available conditions, and no "
            "sustained k-star collapse is detected under the current proxy thresholds."
        )
    else:
        interpretation = (
            "A target-window peak near 3.13 THz is detected, but the spectral audit records a sustained "
            f"recoverability crossing at {k_star_level}. This remains a proxy result, not a materials theory claim."
        )
    return features, interpretation


def run_spectral_feature_audit(
    frequency_spectra: list[dict[str, Any]],
    *,
    time_traces: list[dict[str, Any]] | None = None,
    output_dir: Path,
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
    search_window_THz: float = SEARCH_WINDOW_THZ,
    stft_count: int = 0,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_rows = [
        _audit_frequency_record(
            record,
            target_frequency_THz=target_frequency_THz,
            search_window_THz=search_window_THz,
        )
        for record in frequency_spectra
    ]
    derived_rows = [
        _audit_time_trace_record(
            record,
            target_frequency_THz=target_frequency_THz,
            search_window_THz=search_window_THz,
        )
        for record in (time_traces or [])
    ]
    raw_rows.extend(derived_rows)
    rows = [row for row in raw_rows if row is not None]
    missing_notes: list[str] = []
    plot_notes: list[str] = []

    if stft_count == 0:
        missing_notes.append("No STFT maps parsed in the current data load.")
    derived_count = sum(1 for row in rows if row.get("_derived_from_time_trace"))
    if derived_count:
        missing_notes.append(f"{derived_count} audited spectra are FFT-derived from parsed time traces.")
    if not rows:
        summary = {
            "status": "NO_SPECTRAL_DATA",
            "target_frequency_THz": target_frequency_THz,
            "search_window_THz": search_window_THz,
            "records_audited": 0,
            "peaks_found": 0,
            "standout_features": [],
            "k_star_detected": False,
            "k_star_level": None,
            "visible_failure_detected": False,
            "bounded_interpretation": "No usable spectral records were available for second-pass audit.",
            "missing_data_notes": missing_notes,
            "plots": {},
        }
        _write_json(output_dir / "spectral_feature_summary.json", summary)
        _write_csv(output_dir / "spectral_feature_summary.csv", [], SUMMARY_FIELDNAMES)
        _write_csv(output_dir / "peak_stability_by_condition.csv", [], GROUP_FIELDNAMES)
        return summary

    _apply_group_recoverability(rows)
    distance_decay = audit_distance_decay(rows)
    group_rows = audit_condition_groups(rows)
    k_star_level = _detect_any_k_star(rows)
    standout_features, interpretation = summarize_standout_features(
        rows,
        k_star_level=k_star_level,
        distance_decay=distance_decay,
        target_frequency_THz=target_frequency_THz,
    )

    for row in rows:
        if row.get("linewidth_proxy_THz") is None:
            missing_notes.append(
                f"Linewidth proxy unavailable for {row.get('sample_id')} {row.get('condition_label')}: "
                f"{row.get('_linewidth_note')}"
            )
        if row.get("peak_to_background") is None:
            missing_notes.append(
                f"Peak-to-background unavailable for {row.get('sample_id')} {row.get('condition_label')}: "
                f"{row.get('_background_note')}"
            )

    plot_results = _write_plots(rows, output_dir)
    plot_notes.extend(note for note in plot_results.values() if isinstance(note, str) and note)
    missing_notes.extend(plot_notes)

    peaks_found = sum(1 for row in rows if row.get("peak_found"))
    visible_failure_detected = any(
        (not row.get("peak_found")) or ((_float_or_none(row.get("spectral_recoverability")) or 1.0) < 0.40)
        for row in rows
    )
    status = "PASS"
    if peaks_found == 0:
        status = "NO_SPECTRAL_DATA"
    elif peaks_found < len(rows) or any(row.get("linewidth_proxy_THz") is None for row in rows):
        status = "PARTIAL"

    summary = {
        "status": status,
        "target_frequency_THz": target_frequency_THz,
        "search_window_THz": search_window_THz,
        "records_audited": len(rows),
        "parsed_spectral_records": len(frequency_spectra),
        "time_trace_derived_spectral_records": derived_count,
        "peaks_found": peaks_found,
        "standout_features": standout_features,
        "k_star_detected": k_star_level is not None,
        "k_star_level": k_star_level,
        "visible_failure_detected": visible_failure_detected,
        "bounded_interpretation": interpretation,
        "missing_data_notes": sorted(set(missing_notes)),
        "distance_decay": distance_decay,
        "mean_peak_frequency_THz": _mean_or_none(_finite_values(row.get("peak_frequency_THz") for row in rows)),
        "mean_spectral_recoverability": _mean_or_none(
            _finite_values(row.get("spectral_recoverability") for row in rows)
        ),
        "plots": {
            "spectral_peak_map": str(output_dir / "spectral_peak_map.png")
            if (output_dir / "spectral_peak_map.png").exists()
            else None,
            "peak_frequency_stability": str(output_dir / "peak_frequency_stability.png")
            if (output_dir / "peak_frequency_stability.png").exists()
            else None,
            "peak_amplitude_retention": str(output_dir / "peak_amplitude_retention.png")
            if (output_dir / "peak_amplitude_retention.png").exists()
            else None,
            "linewidth_proxy_by_condition": str(output_dir / "linewidth_proxy_by_condition.png")
            if (output_dir / "linewidth_proxy_by_condition.png").exists()
            else None,
            "spectral_recoverability_audit": str(output_dir / "spectral_recoverability_audit.png")
            if (output_dir / "spectral_recoverability_audit.png").exists()
            else None,
        },
    }

    _write_csv(output_dir / "spectral_feature_summary.csv", rows, SUMMARY_FIELDNAMES)
    _write_csv(output_dir / "peak_stability_by_condition.csv", group_rows, GROUP_FIELDNAMES)
    _write_json(output_dir / "spectral_feature_summary.json", summary)
    return summary


def _audit_frequency_record(
    record: dict[str, Any],
    *,
    target_frequency_THz: float,
    search_window_THz: float,
) -> dict[str, Any] | None:
    source_file, sheet_name, series_label = _parse_source_location(str(record.get("source_file", "")))
    metadata = dict(record.get("metadata") or {})
    distance_um = _float_or_none(record.get("distance_um") if record.get("distance_um") is not None else metadata.get("distance_um"))
    temperature_K = _float_or_none(metadata.get("temperature_K"))
    thickness_nm = _float_or_none(metadata.get("thickness_nm"))
    fluence = _float_or_none(metadata.get("fluence_mJ_cm2"))
    condition_label = _condition_label(record, series_label)
    row: dict[str, Any] = {
        "source_file": source_file,
        "sheet_name": sheet_name,
        "sample_id": str(record.get("sample_id", "unknown")),
        "condition_label": condition_label,
        "distance_um": distance_um,
        "temperature_K": temperature_K,
        "thickness_nm": thickness_nm,
        "fluence_mJ_cm2": fluence,
        "target_frequency_THz": target_frequency_THz,
        "_baseline_group_key": _baseline_group_key(record, sheet_name, series_label),
    }

    if record.get("kind") == "peak_record":
        peak_frequency = _float_or_none(record.get("frequency_THz"))
        peak_amplitude = _float_or_none(record.get("peak_amplitude"))
        found = (
            peak_frequency is not None
            and peak_amplitude is not None
            and abs(peak_frequency - target_frequency_THz) <= search_window_THz
        )
        row.update(
            {
                "peak_found": found,
                "peak_frequency_THz": peak_frequency,
                "peak_amplitude": peak_amplitude,
                "linewidth_proxy_THz": None,
                "peak_to_background": None,
                "frequency_drift_from_target_THz": None
                if peak_frequency is None
                else peak_frequency - target_frequency_THz,
                "_linewidth_note": "peak_record_has_no_spectrum_for_linewidth",
                "_background_note": "peak_record_has_no_spectrum_for_background",
            }
        )
        return row

    frequency = np.asarray(record.get("frequency_THz"), dtype=float)
    amplitude = np.asarray(record.get("amplitude"), dtype=float)
    if frequency.size == 0 or amplitude.size == 0:
        return None
    peak = detect_target_peak(
        frequency,
        amplitude,
        target_frequency_THz=target_frequency_THz,
        search_window_THz=search_window_THz,
    )
    linewidth, linewidth_note = estimate_linewidth_proxy(
        frequency,
        amplitude,
        peak.get("peak_index"),
        target_frequency_THz=target_frequency_THz,
        search_window_THz=search_window_THz,
    )
    peak_to_background, background_note = estimate_peak_to_background(
        frequency,
        amplitude,
        peak.get("peak_amplitude"),
        target_frequency_THz=target_frequency_THz,
        search_window_THz=search_window_THz,
    )
    peak_frequency = _float_or_none(peak.get("peak_frequency_THz"))
    row.update(
        {
            "peak_found": bool(peak.get("peak_found")),
            "peak_frequency_THz": peak_frequency,
            "peak_amplitude": _float_or_none(peak.get("peak_amplitude")),
            "linewidth_proxy_THz": linewidth,
            "peak_to_background": peak_to_background,
            "frequency_drift_from_target_THz": None
            if peak_frequency is None
            else peak_frequency - target_frequency_THz,
            "_linewidth_note": linewidth_note,
            "_background_note": background_note,
        }
    )
    return row


def _audit_time_trace_record(
    record: dict[str, Any],
    *,
    target_frequency_THz: float,
    search_window_THz: float,
) -> dict[str, Any] | None:
    time_ps = np.asarray(record.get("time_ps"), dtype=float)
    signal = np.asarray(record.get("signal"), dtype=float)
    spectrum = _spectrum_from_time_trace(time_ps, signal)
    if spectrum is None:
        return None
    frequency, amplitude = spectrum
    pseudo_record = {
        "kind": "frequency_spectrum",
        "source_file": record.get("source_file", ""),
        "sample_id": record.get("sample_id", "unknown"),
        "condition_type": record.get("condition_type", "dataset"),
        "condition_value": record.get("condition_value", 0.0),
        "distance_um": record.get("distance_um"),
        "metadata": record.get("metadata") or {},
        "frequency_THz": frequency,
        "amplitude": amplitude,
    }
    row = _audit_frequency_record(
        pseudo_record,
        target_frequency_THz=target_frequency_THz,
        search_window_THz=search_window_THz,
    )
    if row is not None:
        row["_derived_from_time_trace"] = True
    return row


def _apply_group_recoverability(rows: list[dict[str, Any]]) -> None:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row.get("_baseline_group_key", "unknown")), []).append(row)

    for group_rows in grouped.values():
        ordered = sorted(group_rows, key=_condition_sort_key)
        baseline = next((row for row in ordered if row.get("peak_found")), ordered[0])
        baseline_amplitude = _float_or_none(baseline.get("peak_amplitude"))
        baseline_linewidth = _float_or_none(baseline.get("linewidth_proxy_THz"))
        baseline_p2b = _float_or_none(baseline.get("peak_to_background"))
        for row in ordered:
            peak_amplitude = _float_or_none(row.get("peak_amplitude"))
            peak_frequency = _float_or_none(row.get("peak_frequency_THz"))
            linewidth = _float_or_none(row.get("linewidth_proxy_THz"))
            peak_to_background = _float_or_none(row.get("peak_to_background"))
            if baseline_amplitude is not None and baseline_amplitude > 1.0e-12 and peak_amplitude is not None:
                amplitude_retention = min(max(peak_amplitude / baseline_amplitude, 0.0), 1.0)
            else:
                amplitude_retention = None
            frequency_stability = None
            if peak_frequency is not None:
                frequency_stability = _clip01(1.0 - abs(peak_frequency - TARGET_FREQUENCY_THZ) / SEARCH_WINDOW_THZ)
            linewidth_stability = None
            if baseline_linewidth is not None and linewidth is not None and baseline_linewidth > 1.0e-12:
                linewidth_stability = _clip01(1.0 - max(linewidth - baseline_linewidth, 0.0) / max(baseline_linewidth, 0.03))
            p2b_stability = None
            if baseline_p2b is not None and peak_to_background is not None and baseline_p2b > 1.0e-12:
                p2b_stability = min(max(peak_to_background / baseline_p2b, 0.0), 1.0)
            row["amplitude_retention"] = amplitude_retention
            row["frequency_stability"] = frequency_stability
            row["linewidth_stability"] = linewidth_stability
            row["_peak_to_background_stability"] = p2b_stability
            row["spectral_recoverability"] = compute_spectral_recoverability(
                amplitude_retention=amplitude_retention,
                frequency_stability=frequency_stability,
                linewidth_stability=linewidth_stability,
                peak_to_background_stability=p2b_stability,
            )


def _detect_any_k_star(rows: list[dict[str, Any]]) -> Any | None:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row.get("_baseline_group_key", "unknown")), []).append(row)
    for group_rows in grouped.values():
        ordered = sorted(group_rows, key=_condition_sort_key)
        values = [_float_or_none(row.get("spectral_recoverability")) for row in ordered]
        labels = [row.get("condition_label") for row in ordered]
        clean = [(value, label) for value, label in zip(values, labels) if value is not None]
        for index in range(len(clean)):
            window = clean[index : index + SUSTAIN_STEPS]
            if len(window) < SUSTAIN_STEPS:
                continue
            if all(value < COLLAPSE_THRESHOLD for value, _ in window):
                return window[0][1]
    return None


def _write_plots(rows: list[dict[str, Any]], output_dir: Path) -> dict[str, str | None]:
    notes: dict[str, str | None] = {}
    notes["spectral_peak_map"] = _plot_peak_map(rows, output_dir / "spectral_peak_map.png")
    notes["peak_frequency_stability"] = _plot_group_metric(
        rows,
        "frequency_stability",
        output_dir / "peak_frequency_stability.png",
        "peak frequency stability",
    )
    notes["peak_amplitude_retention"] = _plot_group_metric(
        rows,
        "amplitude_retention",
        output_dir / "peak_amplitude_retention.png",
        "peak amplitude retention",
    )
    notes["linewidth_proxy_by_condition"] = _plot_group_metric(
        rows,
        "linewidth_proxy_THz",
        output_dir / "linewidth_proxy_by_condition.png",
        "linewidth proxy THz",
    )
    notes["spectral_recoverability_audit"] = _plot_group_metric(
        rows,
        "spectral_recoverability",
        output_dir / "spectral_recoverability_audit.png",
        "spectral recoverability",
    )
    return notes


def _spectrum_from_time_trace(time_ps: np.ndarray, signal: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    time = np.asarray(time_ps, dtype=float)
    values = np.asarray(signal, dtype=float)
    finite = np.isfinite(time) & np.isfinite(values)
    time = time[finite]
    values = values[finite]
    if time.size < 12 or values.size != time.size:
        return None
    order = np.argsort(time)
    time = time[order]
    values = values[order] - float(np.nanmean(values[order]))
    dt = float(np.nanmedian(np.diff(time)))
    if not math.isfinite(dt) or dt <= 0.0:
        return None
    window = np.hanning(values.size)
    frequency = np.fft.rfftfreq(values.size, d=dt)
    amplitude = np.abs(np.fft.rfft(values * window))
    return frequency.astype(float), amplitude.astype(float)


def _plot_peak_map(rows: list[dict[str, Any]], path: Path) -> str | None:
    points = []
    for index, row in enumerate(rows):
        peak_frequency = _float_or_none(row.get("peak_frequency_THz"))
        recoverability = _float_or_none(row.get("spectral_recoverability"))
        amplitude = _float_or_none(row.get("peak_amplitude"))
        if peak_frequency is None or amplitude is None:
            continue
        points.append((index, peak_frequency, amplitude, recoverability if recoverability is not None else 0.0))
    if not points:
        return "spectral_peak_map.png skipped: no detected peak points"
    amplitudes = [point[2] for point in points]
    max_amp = max(amplitudes) if amplitudes else 1.0
    sizes = [35.0 + 115.0 * (amp / max(max_amp, 1.0e-12)) for amp in amplitudes]
    plt.figure(figsize=(7.4, 4.4))
    scatter = plt.scatter(
        [point[0] for point in points],
        [point[1] for point in points],
        s=sizes,
        c=[point[3] for point in points],
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
        edgecolors="black",
        linewidths=0.35,
    )
    plt.axhline(TARGET_FREQUENCY_THZ, color="black", linestyle="--", linewidth=1.0)
    plt.xlabel("audited spectral record")
    plt.ylabel("peak frequency THz")
    plt.colorbar(scatter, label="spectral recoverability")
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()
    return None


def _plot_group_metric(rows: list[dict[str, Any]], metric: str, path: Path, ylabel: str) -> str | None:
    series: dict[str, list[tuple[float, float]]] = {}
    fallback_index = 0
    for row in rows:
        value = _float_or_none(row.get(metric))
        if value is None:
            continue
        x_value = _condition_numeric(row)
        if x_value is None:
            x_value = float(fallback_index)
            fallback_index += 1
        label = str(row.get("_baseline_group_key", row.get("sample_id", "sample")))
        series.setdefault(label, []).append((x_value, value))
    if not series:
        return f"{path.name} skipped: no usable {metric} values"
    plt.figure(figsize=(7.4, 4.4))
    for label, values in sorted(series.items()):
        ordered = sorted(values)
        short_label = label.split("|")[0][-28:]
        plt.plot([x for x, _ in ordered], [y for _, y in ordered], marker="o", linewidth=1.4, label=short_label)
    if metric in {"frequency_stability", "amplitude_retention", "spectral_recoverability"}:
        plt.ylim(0.0, 1.05)
    plt.xlabel("condition value or record index")
    plt.ylabel(ylabel)
    if len(series) <= 7:
        plt.legend(frameon=False, fontsize=7)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()
    return None


def _missing_peak(reason: str) -> dict[str, Any]:
    return {
        "peak_frequency_THz": None,
        "peak_amplitude": None,
        "peak_index": None,
        "peak_found": False,
        "reason": reason,
    }


def _half_max_crossing(
    frequency: np.ndarray,
    adjusted_amplitude: np.ndarray,
    peak_index: int,
    half_height: float,
    *,
    direction: int,
) -> float | None:
    index = peak_index
    while 0 <= index + direction < adjusted_amplitude.size:
        next_index = index + direction
        if adjusted_amplitude[next_index] <= half_height:
            x0 = float(frequency[index])
            x1 = float(frequency[next_index])
            y0 = float(adjusted_amplitude[index])
            y1 = float(adjusted_amplitude[next_index])
            if abs(y1 - y0) <= 1.0e-12:
                return x1
            fraction = (half_height - y0) / (y1 - y0)
            return x0 + fraction * (x1 - x0)
        index = next_index
    return None


def _parse_source_location(source_file: str) -> tuple[str, str | None, str | None]:
    base = source_file
    sheet = None
    series = None
    if "#" in source_file:
        base, remainder = source_file.split("#", 1)
        if ":" in remainder:
            sheet, series = remainder.split(":", 1)
        else:
            sheet = remainder
    return base, sheet, series


def _condition_label(record: dict[str, Any], series_label: str | None) -> str:
    condition_type = str(record.get("condition_type", "dataset"))
    condition_value = record.get("condition_value")
    if condition_type and condition_type != "dataset":
        return f"{condition_type}={condition_value}"
    if series_label:
        return series_label
    return f"{condition_type}={condition_value}"


def _baseline_group_key(record: dict[str, Any], sheet_name: str | None, series_label: str | None) -> str:
    sample_id = str(record.get("sample_id", "unknown"))
    condition_type = str(record.get("condition_type", "dataset"))
    source_file, _, _ = _parse_source_location(str(record.get("source_file", "")))
    if condition_type == "dataset":
        return "|".join([sample_id, source_file, sheet_name or "no_sheet", series_label or "dataset"])
    return "|".join([sample_id, source_file, sheet_name or "no_sheet", condition_type])


def _condition_group_key(row: dict[str, Any]) -> str:
    parts = [
        str(row.get("sample_id", "unknown")),
        str(row.get("source_file", "source")),
        str(row.get("sheet_name", "sheet")),
        str(row.get("condition_label", "condition")),
    ]
    for key in ("distance_um", "temperature_K", "thickness_nm", "fluence_mJ_cm2"):
        value = row.get(key)
        if value not in (None, ""):
            parts.append(f"{key}={value}")
    return "|".join(parts)


def _condition_sort_key(row: dict[str, Any]) -> tuple[int, float | str]:
    numeric = _condition_numeric(row)
    if numeric is not None:
        return (0, numeric)
    return (1, str(row.get("condition_label", "")))


def _condition_numeric(row: dict[str, Any]) -> float | None:
    for key in ("distance_um", "temperature_K", "thickness_nm", "fluence_mJ_cm2"):
        numeric = _float_or_none(row.get(key))
        if numeric is not None:
            return numeric
    label = str(row.get("condition_label", ""))
    if "=" in label:
        numeric = _float_or_none(label.split("=", 1)[1])
        if numeric is not None:
            return numeric
    return None


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _csv_value(row.get(field)) for field in fieldnames})


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


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


def _finite_values(values: Any) -> list[float]:
    finite: list[float] = []
    for value in values:
        numeric = _float_or_none(value)
        if numeric is not None:
            finite.append(numeric)
    return finite


def _mean_or_none(values: list[float]) -> float | None:
    if not values:
        return None
    return float(sum(values) / len(values))


def _std_or_none(values: list[float]) -> float | None:
    if not values:
        return None
    mean = sum(values) / len(values)
    return float(math.sqrt(sum((value - mean) ** 2 for value in values) / len(values)))
