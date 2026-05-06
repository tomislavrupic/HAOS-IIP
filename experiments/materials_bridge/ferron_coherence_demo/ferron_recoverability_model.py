#!/usr/bin/env python3
"""HAOS-style recoverability proxies for ferron/NbOI2 sidecar data.

The metrics here are not materials theory. They are bounded telemetry proxies
for asking whether a narrow-band coherent readout remains recoverable under
distance, delay, or other perturbation axes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np


TARGET_FREQUENCY_THZ = 3.13
TARGET_WINDOW_THZ = 0.38
COLLAPSE_THRESHOLD = 0.70
VISIBLE_AMPLITUDE_THRESHOLD = 0.40
VISIBLE_LINEWIDTH_RATIO_THRESHOLD = 2.50
SUSTAIN_STEPS = 2


@dataclass(frozen=True)
class Measurement:
    source_mode: str
    source_file: str
    sample_id: str
    condition_type: str
    condition_value: Any
    distance_um: float | None
    time_ps: float | None
    frequency_THz: float | None
    peak_amplitude: float
    peak_frequency_THz: float | None
    linewidth_proxy: float | None
    coherence_proxy: float | None
    peak_to_background_proxy: float | None


def normalize_signal(signal: np.ndarray) -> np.ndarray:
    values = np.asarray(signal, dtype=float)
    if values.size == 0:
        return values
    centered = values - float(np.nanmean(values))
    scale = float(np.nanmax(np.abs(centered))) if centered.size else 0.0
    if not math.isfinite(scale) or scale <= 1.0e-12:
        return np.zeros_like(centered)
    return centered / scale


def estimate_peak_amplitude(
    frequency_THz: np.ndarray,
    amplitude: np.ndarray,
    *,
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
    window_THz: float = TARGET_WINDOW_THZ,
) -> float:
    frequency = np.asarray(frequency_THz, dtype=float)
    values = np.abs(np.asarray(amplitude, dtype=float))
    if frequency.size == 0 or values.size == 0 or frequency.size != values.size:
        return float("nan")
    mask = np.abs(frequency - target_frequency_THz) <= window_THz
    if not np.any(mask):
        return float("nan")
    background = _background_level(frequency, values, target_frequency_THz, window_THz)
    return max(float(np.nanmax(values[mask]) - background), 0.0)


def estimate_peak_frequency(
    frequency_THz: np.ndarray,
    amplitude: np.ndarray,
    *,
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
    window_THz: float = TARGET_WINDOW_THZ,
) -> float:
    frequency = np.asarray(frequency_THz, dtype=float)
    values = np.abs(np.asarray(amplitude, dtype=float))
    if frequency.size == 0 or values.size == 0 or frequency.size != values.size:
        return float("nan")
    mask = np.abs(frequency - target_frequency_THz) <= window_THz
    if not np.any(mask):
        return float("nan")
    local_indices = np.where(mask)[0]
    peak_index = local_indices[int(np.nanargmax(values[mask]))]
    return float(frequency[peak_index])


def estimate_linewidth_proxy(
    frequency_THz: np.ndarray,
    amplitude: np.ndarray,
    *,
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
    window_THz: float = TARGET_WINDOW_THZ,
) -> float:
    frequency = np.asarray(frequency_THz, dtype=float)
    values = np.abs(np.asarray(amplitude, dtype=float))
    if frequency.size < 4 or values.size != frequency.size:
        return float("nan")
    mask = np.abs(frequency - target_frequency_THz) <= window_THz
    if not np.any(mask):
        return float("nan")
    local_frequency = frequency[mask]
    local_values = values[mask]
    background = float(np.nanmin(local_values))
    adjusted = np.maximum(local_values - background, 0.0)
    peak = float(np.nanmax(adjusted))
    if peak <= 1.0e-12:
        return float("nan")
    half = 0.5 * peak
    peak_index = int(np.nanargmax(adjusted))
    above = adjusted >= half
    left = peak_index
    while left > 0 and above[left - 1]:
        left -= 1
    right = peak_index
    while right < above.size - 1 and above[right + 1]:
        right += 1
    width = float(local_frequency[right] - local_frequency[left])
    if width <= 0.0 and local_frequency.size > 1:
        width = float(np.nanmedian(np.diff(np.sort(local_frequency))))
    return max(width, 0.0)


def estimate_coherence_proxy(
    time_ps: np.ndarray,
    signal: np.ndarray,
    *,
    target_frequency_THz: float = TARGET_FREQUENCY_THZ,
) -> float:
    time = np.asarray(time_ps, dtype=float)
    values = normalize_signal(np.asarray(signal, dtype=float))
    if time.size < 12 or values.size != time.size:
        return float("nan")
    spectrum_frequency, spectrum_amplitude = _spectrum_from_trace(time, values)
    if spectrum_frequency.size == 0:
        return float("nan")
    peak = estimate_peak_amplitude(spectrum_frequency, spectrum_amplitude, target_frequency_THz=target_frequency_THz)
    background = _background_level(spectrum_frequency, np.abs(spectrum_amplitude), target_frequency_THz, TARGET_WINDOW_THZ)
    spectral_clarity = peak / max(peak + background, 1.0e-12)

    midpoint = values.size // 2
    early = float(np.sqrt(np.mean(values[:midpoint] ** 2)))
    late = float(np.sqrt(np.mean(values[midpoint:] ** 2)))
    envelope_survival = late / max(early, 1.0e-12)
    return _clip01(0.65 * spectral_clarity + 0.35 * min(envelope_survival, 1.0))


def compute_recoverability(
    peak_amplitude_normalized: float,
    peak_frequency_stability: float,
    linewidth_stability: float,
) -> float:
    if not math.isfinite(peak_amplitude_normalized):
        return 0.0
    return _clip01(peak_amplitude_normalized) * _clip01(peak_frequency_stability) * _clip01(linewidth_stability)


def compute_delta_persistence(recoverability_values: list[float]) -> list[float]:
    deltas: list[float] = []
    previous: float | None = None
    for value in recoverability_values:
        if previous is None:
            deltas.append(0.0)
        else:
            deltas.append(float(value - previous))
        previous = value
    return deltas


def detect_k_star(
    recoverability_values: list[float],
    condition_values: list[Any],
    *,
    collapse_threshold: float = COLLAPSE_THRESHOLD,
    sustain_steps: int = SUSTAIN_STEPS,
) -> Any | None:
    if sustain_steps <= 0:
        sustain_steps = 1
    for index in range(len(recoverability_values)):
        window = recoverability_values[index : index + sustain_steps]
        if len(window) < sustain_steps:
            continue
        if all(value < collapse_threshold for value in window):
            return condition_values[index]
    return None


def compute_safety_margin(current_condition: Any, k_star: Any | None, all_conditions: list[Any]) -> float | None:
    if k_star is None:
        return None
    current_numeric = _float_or_none(current_condition)
    k_numeric = _float_or_none(k_star)
    if current_numeric is not None and k_numeric is not None:
        return float(k_numeric - current_numeric)
    try:
        return float(all_conditions.index(k_star) - all_conditions.index(current_condition))
    except ValueError:
        return None


def compute_confidence(
    delta_persistence_so_far: list[float],
    *,
    peak_to_background_proxy: float | None = None,
) -> float:
    if peak_to_background_proxy is not None and math.isfinite(peak_to_background_proxy):
        clarity = _clip01(peak_to_background_proxy)
    else:
        clarity = 0.0
    if delta_persistence_so_far:
        strongest_drop = abs(min(delta_persistence_so_far))
    else:
        strongest_drop = 0.0
    return _clip01(max(clarity, strongest_drop))


def detect_visible_failure(
    *,
    peak_amplitude_normalized: float,
    linewidth_proxy: float | None,
    baseline_linewidth_proxy: float | None,
    peak_frequency_THz: float | None,
) -> bool:
    if peak_frequency_THz is None or not math.isfinite(peak_frequency_THz):
        return True
    if not math.isfinite(peak_amplitude_normalized) or peak_amplitude_normalized < VISIBLE_AMPLITUDE_THRESHOLD:
        return True
    if (
        linewidth_proxy is not None
        and baseline_linewidth_proxy is not None
        and math.isfinite(linewidth_proxy)
        and math.isfinite(baseline_linewidth_proxy)
        and baseline_linewidth_proxy > 1.0e-12
        and linewidth_proxy / baseline_linewidth_proxy > VISIBLE_LINEWIDTH_RATIO_THRESHOLD
    ):
        return True
    return False


def run_ferron_sweep(
    *,
    time_traces: list[dict[str, Any]] | None = None,
    frequency_spectra: list[dict[str, Any]] | None = None,
    stft_amplitudes: list[dict[str, Any]] | None = None,
    source_mode: str = "REAL_DATA_LOADED",
) -> list[dict[str, Any]]:
    measurements: list[Measurement] = []
    for record in frequency_spectra or []:
        measurements.extend(_measure_from_frequency_record(record, source_mode))
    for record in time_traces or []:
        measurement = _measure_from_time_trace(record, source_mode)
        if measurement is not None:
            measurements.append(measurement)
    for record in stft_amplitudes or []:
        measurements.extend(_measure_from_stft_record(record, source_mode))

    grouped: dict[tuple[str, str, str, str], list[Measurement]] = {}
    for measurement in measurements:
        key = (
            measurement.source_mode,
            measurement.sample_id,
            measurement.condition_type,
            _source_family(measurement.source_file),
        )
        grouped.setdefault(key, []).append(measurement)

    rows: list[dict[str, Any]] = []
    for group_measurements in grouped.values():
        rows.extend(_compute_group_rows(group_measurements))
    return rows


def _measure_from_frequency_record(record: dict[str, Any], source_mode: str) -> list[Measurement]:
    if record.get("kind") == "peak_record":
        peak = _float_or_nan(record.get("peak_amplitude"))
        frequency = _float_or_none(record.get("frequency_THz"))
        coherence_time = _float_or_none(record.get("coherence_time_ps"))
        linewidth = (1.0 / coherence_time) if coherence_time and coherence_time > 0 else None
        return [
            Measurement(
                source_mode=source_mode,
                source_file=str(record.get("source_file", "")),
                sample_id=str(record.get("sample_id", "unknown")),
                condition_type=str(record.get("condition_type", "condition_index")),
                condition_value=record.get("condition_value", 0.0),
                distance_um=_float_or_none(record.get("distance_um")),
                time_ps=_time_from_condition(record),
                frequency_THz=frequency,
                peak_amplitude=peak,
                peak_frequency_THz=frequency,
                linewidth_proxy=linewidth,
                coherence_proxy=None,
                peak_to_background_proxy=None,
            )
        ]

    frequency_axis = np.asarray(record.get("frequency_THz"), dtype=float)
    amplitude = np.asarray(record.get("amplitude"), dtype=float)
    peak = estimate_peak_amplitude(frequency_axis, amplitude)
    peak_frequency = estimate_peak_frequency(frequency_axis, amplitude)
    linewidth = estimate_linewidth_proxy(frequency_axis, amplitude)
    background = _background_level(frequency_axis, np.abs(amplitude), TARGET_FREQUENCY_THZ, TARGET_WINDOW_THZ)
    clarity = peak / max(peak + background, 1.0e-12) if math.isfinite(peak) else 0.0
    return [
        Measurement(
            source_mode=source_mode,
            source_file=str(record.get("source_file", "")),
            sample_id=str(record.get("sample_id", "unknown")),
            condition_type=str(record.get("condition_type", "dataset")),
            condition_value=record.get("condition_value", 0.0),
            distance_um=_float_or_none(record.get("distance_um")),
            time_ps=_time_from_condition(record),
            frequency_THz=None,
            peak_amplitude=peak,
            peak_frequency_THz=peak_frequency if math.isfinite(peak_frequency) else None,
            linewidth_proxy=linewidth if math.isfinite(linewidth) else None,
            coherence_proxy=None,
            peak_to_background_proxy=clarity,
        )
    ]


def _measure_from_time_trace(record: dict[str, Any], source_mode: str) -> Measurement | None:
    time = np.asarray(record.get("time_ps"), dtype=float)
    signal = np.asarray(record.get("signal"), dtype=float)
    if time.size < 12 or signal.size != time.size:
        return None
    frequency, amplitude = _spectrum_from_trace(time, signal)
    if frequency.size == 0:
        return None
    peak = estimate_peak_amplitude(frequency, amplitude)
    peak_frequency = estimate_peak_frequency(frequency, amplitude)
    linewidth = estimate_linewidth_proxy(frequency, amplitude)
    coherence = estimate_coherence_proxy(time, signal)
    background = _background_level(frequency, np.abs(amplitude), TARGET_FREQUENCY_THZ, TARGET_WINDOW_THZ)
    clarity = peak / max(peak + background, 1.0e-12) if math.isfinite(peak) else 0.0
    return Measurement(
        source_mode=source_mode,
        source_file=str(record.get("source_file", "")),
        sample_id=str(record.get("sample_id", "unknown")),
        condition_type=str(record.get("condition_type", "dataset")),
        condition_value=record.get("condition_value", 0.0),
        distance_um=_float_or_none(record.get("distance_um")),
        time_ps=_time_from_condition(record),
        frequency_THz=None,
        peak_amplitude=peak,
        peak_frequency_THz=peak_frequency if math.isfinite(peak_frequency) else None,
        linewidth_proxy=linewidth if math.isfinite(linewidth) else None,
        coherence_proxy=coherence if math.isfinite(coherence) else None,
        peak_to_background_proxy=clarity,
    )


def _measure_from_stft_record(record: dict[str, Any], source_mode: str) -> list[Measurement]:
    time = np.asarray(record.get("time_ps"), dtype=float)
    frequency = np.asarray(record.get("frequency_THz"), dtype=float)
    amplitude = np.asarray(record.get("amplitude"), dtype=float)
    if time.size < 8 or frequency.size != time.size or amplitude.size != time.size:
        return []

    measurements: list[Measurement] = []
    for time_value in sorted(set(float(value) for value in time)):
        mask = np.isclose(time, time_value)
        if int(np.sum(mask)) < 4:
            continue
        freq_slice = frequency[mask]
        amp_slice = amplitude[mask]
        peak = estimate_peak_amplitude(freq_slice, amp_slice)
        peak_frequency = estimate_peak_frequency(freq_slice, amp_slice)
        linewidth = estimate_linewidth_proxy(freq_slice, amp_slice)
        background = _background_level(freq_slice, np.abs(amp_slice), TARGET_FREQUENCY_THZ, TARGET_WINDOW_THZ)
        clarity = peak / max(peak + background, 1.0e-12) if math.isfinite(peak) else 0.0
        measurements.append(
            Measurement(
                source_mode=source_mode,
                source_file=str(record.get("source_file", "")),
                sample_id=str(record.get("sample_id", "unknown")),
                condition_type="time_ps",
                condition_value=float(time_value),
                distance_um=_float_or_none(record.get("distance_um")),
                time_ps=float(time_value),
                frequency_THz=None,
                peak_amplitude=peak,
                peak_frequency_THz=peak_frequency if math.isfinite(peak_frequency) else None,
                linewidth_proxy=linewidth if math.isfinite(linewidth) else None,
                coherence_proxy=None,
                peak_to_background_proxy=clarity,
            )
        )
    return measurements


def _compute_group_rows(measurements: list[Measurement]) -> list[dict[str, Any]]:
    ordered = sorted(enumerate(measurements), key=lambda item: _condition_sort_key(item[1].condition_value, item[0]))
    ordered_measurements = [measurement for _, measurement in ordered]
    if not ordered_measurements:
        return []

    baseline = ordered_measurements[0]
    baseline_amp = baseline.peak_amplitude if math.isfinite(baseline.peak_amplitude) else 0.0
    baseline_freq = baseline.peak_frequency_THz or TARGET_FREQUENCY_THZ
    baseline_linewidth = baseline.linewidth_proxy
    baseline_coherence = baseline.coherence_proxy

    prelim: list[dict[str, Any]] = []
    combined_values: list[float] = []
    condition_values: list[Any] = []
    for condition_index, measurement in enumerate(ordered_measurements):
        amp_norm = measurement.peak_amplitude / max(baseline_amp, 1.0e-12)
        amp_norm = max(amp_norm, 0.0) if math.isfinite(amp_norm) else 0.0
        freq_stability = _frequency_stability(measurement.peak_frequency_THz, baseline_freq)
        linewidth_stability = _linewidth_stability(measurement.linewidth_proxy, baseline_linewidth)
        peak_recoverability = compute_recoverability(amp_norm, freq_stability, linewidth_stability)
        propagation = _clip01(amp_norm) if measurement.condition_type == "distance_um" else None
        temporal = None
        if measurement.condition_type in {"time_ps", "delay_ps"}:
            temporal = _clip01(amp_norm)
        elif measurement.coherence_proxy is not None and baseline_coherence is not None and baseline_coherence > 1.0e-12:
            temporal = _clip01(measurement.coherence_proxy / baseline_coherence)

        combined_parts = [peak_recoverability]
        if propagation is not None:
            combined_parts.append(propagation)
        if temporal is not None:
            combined_parts.append(temporal)
        combined = float(np.mean(combined_parts))

        visible_failure = detect_visible_failure(
            peak_amplitude_normalized=amp_norm,
            linewidth_proxy=measurement.linewidth_proxy,
            baseline_linewidth_proxy=baseline_linewidth,
            peak_frequency_THz=measurement.peak_frequency_THz,
        )
        condition_values.append(measurement.condition_value)
        combined_values.append(combined)
        prelim.append(
            {
                "source_mode": measurement.source_mode,
                "sample_id": measurement.sample_id,
                "condition_index": condition_index,
                "condition_type": measurement.condition_type,
                "condition_value": measurement.condition_value,
                "distance_um": measurement.distance_um,
                "time_ps": measurement.time_ps,
                "frequency_THz": measurement.frequency_THz,
                "peak_amplitude": measurement.peak_amplitude,
                "peak_frequency_THz": measurement.peak_frequency_THz,
                "linewidth_proxy": measurement.linewidth_proxy,
                "ferron_peak_recoverability": peak_recoverability,
                "propagation_recoverability": propagation,
                "temporal_recoverability": temporal,
                "combined_recoverability": combined,
                "visible_failure": visible_failure,
                "_peak_to_background_proxy": measurement.peak_to_background_proxy,
            }
        )

    deltas = compute_delta_persistence(combined_values)
    k_star = detect_k_star(combined_values, condition_values)
    for index, row in enumerate(prelim):
        row["delta_persistence"] = deltas[index]
        row["k_star"] = k_star
        row["safety_margin"] = compute_safety_margin(row["condition_value"], k_star, condition_values)
        row["confidence"] = compute_confidence(
            deltas[: index + 1],
            peak_to_background_proxy=row.pop("_peak_to_background_proxy"),
        )
    return prelim


def _spectrum_from_trace(time_ps: np.ndarray, signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    time = np.asarray(time_ps, dtype=float)
    values = normalize_signal(np.asarray(signal, dtype=float))
    if time.size < 8 or values.size != time.size:
        return np.asarray([]), np.asarray([])
    order = np.argsort(time)
    time = time[order]
    values = values[order]
    dt = float(np.nanmedian(np.diff(time)))
    if not math.isfinite(dt) or dt <= 0.0:
        return np.asarray([]), np.asarray([])
    window = np.hanning(values.size)
    spectrum = np.abs(np.fft.rfft(values * window))
    frequency = np.fft.rfftfreq(values.size, d=dt)
    return frequency.astype(float), spectrum.astype(float)


def _background_level(
    frequency_THz: np.ndarray,
    amplitude: np.ndarray,
    target_frequency_THz: float,
    window_THz: float,
) -> float:
    frequency = np.asarray(frequency_THz, dtype=float)
    values = np.abs(np.asarray(amplitude, dtype=float))
    if frequency.size == 0 or values.size != frequency.size:
        return 0.0
    sideband = np.abs(frequency - target_frequency_THz) > window_THz
    if np.any(sideband):
        return float(np.nanmedian(values[sideband]))
    return float(np.nanmin(values))


def _frequency_stability(peak_frequency: float | None, baseline_frequency: float) -> float:
    if peak_frequency is None or not math.isfinite(peak_frequency):
        return 0.0
    scale = max(0.10 * abs(baseline_frequency), 0.05)
    drift = abs(float(peak_frequency) - float(baseline_frequency))
    return _clip01(1.0 - drift / scale)


def _linewidth_stability(linewidth: float | None, baseline_linewidth: float | None) -> float:
    if linewidth is None or baseline_linewidth is None:
        return 1.0
    if not math.isfinite(linewidth) or not math.isfinite(baseline_linewidth):
        return 1.0
    scale = max(baseline_linewidth, 0.03)
    broadening = max(float(linewidth - baseline_linewidth), 0.0)
    return _clip01(1.0 - broadening / scale)


def _clip01(value: float) -> float:
    if not math.isfinite(value):
        return 0.0
    return min(max(float(value), 0.0), 1.0)


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


def _float_or_nan(value: Any) -> float:
    numeric = _float_or_none(value)
    return numeric if numeric is not None else float("nan")


def _time_from_condition(record: dict[str, Any]) -> float | None:
    if record.get("condition_type") in {"time_ps", "delay_ps"}:
        return _float_or_none(record.get("condition_value"))
    return None


def _condition_sort_key(value: Any, fallback_index: int) -> tuple[int, float | str]:
    numeric = _float_or_none(value)
    if numeric is not None:
        return (0, numeric)
    if value is None:
        return (1, float(fallback_index))
    return (2, str(value))


def _source_family(source_file: str) -> str:
    lower = source_file.lower()
    if "#" in lower and ":" in lower:
        lower = lower.split(":", 1)[0]
    if "time" in lower or "trace" in lower:
        return "time_trace"
    if "stft" in lower or "timefreq" in lower:
        return "stft"
    if "spectrum" in lower or "freq" in lower or "fft" in lower:
        return "spectrum"
    return lower
