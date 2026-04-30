#!/usr/bin/env python3
"""Phase 64 GW memory entry gate.

This sidecar is a claim-gated bridge from toy celestial probes toward
strain-like time-series data. It asks a deliberately narrow question:

Can the existing HAOS-style memory/composition telemetry distinguish a
structured strain-derived proxy event from controls?

The answer is not a gravitational-memory detection. The runner labels every
output as a strain-derived proxy and keeps real BMS, soft theorem, celestial
holography, and detector-memory claims OPEN.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from gw_data_loader import StrainRecord, load_strain_record, prepare_record, write_record_metadata, write_record_preview


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(__file__).resolve().parent

CSV_PATH = ROOT / "gw_memory_diagnostics.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "gw_memory_entry_report.md"
FIGURE_DIR = ROOT / "figures"
PREVIEW_PATH = ROOT / "prepared_strain_preview.csv"
METADATA_PATH = ROOT / "prepared_strain_metadata.json"

FIELDNAMES = [
    "sample",
    "sample_role",
    "control_type",
    "strain_proxy_score",
    "memory_proxy_recoverability",
    "event_localization_score",
    "chirp_order_score",
    "composition_alignment_score",
    "permanence_score",
    "baseline_quiet_score",
    "generic_band_energy_score",
    "event_energy_fraction",
    "peak_time_error",
    "centroid_slope_hz",
    "step_profile_correlation",
    "plateau_stability",
]


@dataclass(frozen=True)
class ProxyConfig:
    duration: float = 2.0
    analysis_sample_rate: float = 1024.0
    event_time: float | None = None
    seed: int = 6401
    signal_to_noise: float = 14.0
    f_min: float = 24.0
    f_max: float = 360.0
    event_half_width: float = 0.22
    event_width: float = 0.090
    smoothing_seconds: float = 0.024
    pass_min_proxy_score: float = 0.76
    pass_min_memory_recoverability: float = 0.72
    pass_min_event_localization: float = 0.58
    pass_min_chirp_order: float = 0.56
    pass_min_composition_alignment: float = 0.66
    pass_min_margin: float = 0.14
    pass_min_separation_z: float = 1.50
    pass_max_control_false_positive_rate: float = 0.0


@dataclass(frozen=True)
class ProxyProducts:
    filtered: np.ndarray
    envelope: np.ndarray
    flux: np.ndarray
    cumulative_memory_proxy: np.ndarray
    ideal_gate: np.ndarray


@dataclass(frozen=True)
class SampleSpec:
    sample: str
    sample_role: str
    control_type: str
    strain: np.ndarray


@dataclass(frozen=True)
class SampleSummary:
    sample: str
    sample_role: str
    control_type: str
    strain_proxy_score: float
    memory_proxy_recoverability: float
    event_localization_score: float
    chirp_order_score: float
    composition_alignment_score: float
    permanence_score: float
    baseline_quiet_score: float
    generic_band_energy_score: float
    event_energy_fraction: float
    peak_time_error: float
    centroid_slope_hz: float
    step_profile_correlation: float
    plateau_stability: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 64 GW memory entry-gate proxy.")
    parser.add_argument("--input-file", type=Path, default=None, help="Optional local GWOSC HDF5, CSV/TXT, or NPY strain file.")
    parser.add_argument("--duration", type=float, default=2.0, help="Analysis-window duration in seconds.")
    parser.add_argument("--sample-rate", type=float, default=1024.0, help="Synthetic or text/NPY sample rate.")
    parser.add_argument("--analysis-sample-rate", type=float, default=1024.0, help="Resample analysis segment to this rate.")
    parser.add_argument("--event-time", type=float, default=None, help="Event time in file coordinates; defaults to surrogate event or file midpoint.")
    parser.add_argument("--seed", type=int, default=6401)
    parser.add_argument("--signal-to-noise", type=float, default=14.0, help="Synthetic surrogate analysis-scale signal-to-noise.")
    parser.add_argument("--detector", type=str, default="SYNTH")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated CSV, JSON, markdown report, preview data, and figures.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> ProxyConfig:
    return ProxyConfig(
        duration=float(args.duration),
        analysis_sample_rate=float(args.analysis_sample_rate),
        event_time=None if args.event_time is None else float(args.event_time),
        seed=int(args.seed),
        signal_to_noise=float(args.signal_to_noise),
    )


def normalize_signal(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    arr = arr - float(np.median(arr))
    scale = float(np.std(arr))
    if scale <= 1.0e-12:
        return np.zeros_like(arr)
    return arr / scale


def cosine_similarity(left: np.ndarray, right: np.ndarray) -> float:
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denom <= 1.0e-12:
        return 0.0
    return float(np.dot(left, right) / denom)


def pearson_correlation(left: np.ndarray, right: np.ndarray) -> float:
    return cosine_similarity(left - float(np.mean(left)), right - float(np.mean(right)))


def sigmoid_gate(time: np.ndarray, event_time: float, width: float) -> np.ndarray:
    return 0.5 * (1.0 + np.tanh((time - event_time) / max(width, 1.0e-9)))


def gaussian_gate(time: np.ndarray, center: float, width: float) -> np.ndarray:
    return np.exp(-0.5 * ((time - center) / max(width, 1.0e-9)) ** 2)


def fft_bandpass(strain: np.ndarray, sample_rate: float, f_min: float, f_max: float) -> np.ndarray:
    centered = normalize_signal(strain)
    spectrum = np.fft.rfft(centered)
    freqs = np.fft.rfftfreq(centered.size, d=1.0 / float(sample_rate))
    mask = (freqs >= f_min) & (freqs <= f_max)
    spectrum = spectrum * mask
    filtered = np.fft.irfft(spectrum, n=centered.size)
    return normalize_signal(filtered)


def moving_average(values: np.ndarray, window_count: int) -> np.ndarray:
    window_count = max(1, int(window_count))
    kernel = np.ones(window_count, dtype=float) / float(window_count)
    return np.convolve(values, kernel, mode="same")


def moving_rms(values: np.ndarray, sample_rate: float, seconds: float) -> np.ndarray:
    window = max(3, int(round(seconds * sample_rate)))
    return np.sqrt(np.maximum(moving_average(values * values, window), 0.0))


def event_masks(time: np.ndarray, event_time: float, config: ProxyConfig) -> dict[str, np.ndarray]:
    event_half = config.event_half_width
    masks = {
        "pre": (time >= event_time - 0.78) & (time <= event_time - 0.34),
        "event": (time >= event_time - event_half) & (time <= event_time + event_half * 0.58),
        "post": (time >= event_time + 0.30) & (time <= event_time + 0.72),
        "early": (time >= event_time - 0.20) & (time < event_time - 0.085),
        "middle": (time >= event_time - 0.085) & (time < event_time + 0.020),
        "late": (time >= event_time + 0.020) & (time <= event_time + 0.135),
    }
    for name, mask in list(masks.items()):
        if np.count_nonzero(mask) < 4:
            masks[name] = _nearest_mask(time, event_time, width=0.08 if name != "post" else 0.12)
    return masks


def _nearest_mask(time: np.ndarray, center: float, width: float) -> np.ndarray:
    return np.abs(time - center) <= width


def compute_products(record: StrainRecord, config: ProxyConfig, strain: np.ndarray) -> ProxyProducts:
    time = record.time
    filtered = fft_bandpass(strain, record.sample_rate, config.f_min, config.f_max)
    envelope = moving_rms(filtered, record.sample_rate, config.smoothing_seconds)
    masks = event_masks(time, record.event_time, config)
    baseline = float(np.median(envelope[masks["pre"]])) if np.any(masks["pre"]) else float(np.median(envelope))
    excess = np.maximum(envelope - baseline, 0.0)
    flux = excess * excess
    total = float(np.sum(flux))
    if total <= 1.0e-12:
        cumulative = np.zeros_like(flux)
    else:
        cumulative = np.cumsum(flux) / total
    ideal_gate = sigmoid_gate(time, record.event_time, config.event_width)
    return ProxyProducts(
        filtered=filtered,
        envelope=envelope,
        flux=flux,
        cumulative_memory_proxy=cumulative,
        ideal_gate=ideal_gate,
    )


def spectral_centroid(segment: np.ndarray, sample_rate: float, config: ProxyConfig) -> float:
    segment = normalize_signal(segment)
    if segment.size < 8:
        return 0.0
    window = np.hanning(segment.size)
    power = np.abs(np.fft.rfft(segment * window)) ** 2
    freqs = np.fft.rfftfreq(segment.size, d=1.0 / float(sample_rate))
    mask = (freqs >= config.f_min) & (freqs <= config.f_max)
    power = power[mask]
    freqs = freqs[mask]
    denom = float(np.sum(power))
    if denom <= 1.0e-12:
        return 0.0
    return float(np.sum(freqs * power) / denom)


def chirp_order_score(record: StrainRecord, config: ProxyConfig, filtered: np.ndarray) -> tuple[float, float]:
    masks = event_masks(record.time, record.event_time, config)
    centroids = np.array(
        [
            spectral_centroid(filtered[masks["early"]], record.sample_rate, config),
            spectral_centroid(filtered[masks["middle"]], record.sample_rate, config),
            spectral_centroid(filtered[masks["late"]], record.sample_rate, config),
        ],
        dtype=float,
    )
    diffs = np.diff(centroids)
    monotone = float(np.mean(diffs > 4.0))
    slope = float(centroids[-1] - centroids[0])
    slope_score = float(np.clip(slope / 95.0, 0.0, 1.0))
    corr = pearson_correlation(centroids, np.arange(centroids.size, dtype=float))
    corr_score = float(np.clip(max(corr, 0.0), 0.0, 1.0))
    return float(0.34 * monotone + 0.34 * slope_score + 0.32 * corr_score), slope


def memory_proxy_recoverability(
    record: StrainRecord,
    config: ProxyConfig,
    products: ProxyProducts,
) -> tuple[float, float, float, float, float]:
    masks = event_masks(record.time, record.event_time, config)
    cumulative = products.cumulative_memory_proxy
    corr = pearson_correlation(cumulative, products.ideal_gate)
    corr_score = float(np.clip(max(corr, 0.0), 0.0, 1.0))
    pre_mean = float(np.mean(cumulative[masks["pre"]])) if np.any(masks["pre"]) else float(cumulative[0])
    post_mean = float(np.mean(cumulative[masks["post"]])) if np.any(masks["post"]) else float(cumulative[-1])
    contrast = float(np.clip((post_mean - pre_mean) / 0.82, 0.0, 1.0))
    plateau_std = float(np.std(cumulative[masks["post"]])) if np.any(masks["post"]) else float(np.std(cumulative[-10:]))
    plateau = float(np.clip(1.0 - plateau_std / 0.070, 0.0, 1.0))
    baseline_std = float(np.std(cumulative[masks["pre"]])) if np.any(masks["pre"]) else float(np.std(cumulative[:10]))
    baseline_quiet = float(np.clip(1.0 - baseline_std / 0.045, 0.0, 1.0))
    score = float(0.44 * corr_score + 0.26 * contrast + 0.18 * plateau + 0.12 * baseline_quiet)
    return score, corr, plateau, baseline_quiet, contrast


def event_localization_score(
    record: StrainRecord,
    config: ProxyConfig,
    products: ProxyProducts,
) -> tuple[float, float, float]:
    masks = event_masks(record.time, record.event_time, config)
    total = max(float(np.sum(products.flux)), 1.0e-12)
    event_fraction = float(np.sum(products.flux[masks["event"]]) / total)
    fraction_score = float(np.clip((event_fraction - 0.22) / 0.52, 0.0, 1.0))
    peak_time = float(record.time[int(np.argmax(products.envelope))])
    peak_error = abs(peak_time - float(record.event_time))
    peak_score = float(np.clip(1.0 - peak_error / 0.22, 0.0, 1.0))
    return float(0.68 * fraction_score + 0.32 * peak_score), event_fraction, peak_error


def composition_alignment_score(record: StrainRecord, config: ProxyConfig, products: ProxyProducts) -> float:
    derivative = np.gradient(products.cumulative_memory_proxy, record.time)
    event_shape = gaussian_gate(record.time, record.event_time - 0.035, config.event_half_width * 0.50)
    event_shape = event_shape / max(float(np.sum(event_shape)), 1.0e-12)
    derivative = derivative / max(float(np.sum(np.abs(derivative))), 1.0e-12)
    corr = pearson_correlation(derivative, event_shape)
    corr_score = float(np.clip(max(corr, 0.0), 0.0, 1.0))
    derivative_peak_time = float(record.time[int(np.argmax(derivative))])
    peak_score = float(np.clip(1.0 - abs(derivative_peak_time - record.event_time) / 0.20, 0.0, 1.0))
    return float(0.72 * corr_score + 0.28 * peak_score)


def generic_band_energy_score(record: StrainRecord, config: ProxyConfig, filtered: np.ndarray) -> float:
    raw_power = np.abs(np.fft.rfft(normalize_signal(record.strain))) ** 2
    filtered_power = np.abs(np.fft.rfft(filtered)) ** 2
    denom = max(float(np.sum(raw_power)), 1.0e-12)
    fraction = float(np.sum(filtered_power) / denom)
    return float(np.clip(fraction / 0.70, 0.0, 1.0))


def evaluate_sample(record: StrainRecord, config: ProxyConfig, sample: SampleSpec) -> SampleSummary:
    products = compute_products(record, config, sample.strain)
    memory_score, step_corr, plateau, baseline_quiet, _contrast = memory_proxy_recoverability(record, config, products)
    event_score, event_fraction, peak_error = event_localization_score(record, config, products)
    chirp_score, centroid_slope = chirp_order_score(record, config, products.filtered)
    composition = composition_alignment_score(record, config, products)
    generic = generic_band_energy_score(record, config, products.filtered)
    permanence = float(0.70 * plateau + 0.30 * baseline_quiet)
    proxy_score = float(
        0.31 * memory_score
        + 0.25 * chirp_score
        + 0.22 * event_score
        + 0.17 * composition
        + 0.03 * permanence
        + 0.02 * generic
    )
    return SampleSummary(
        sample=sample.sample,
        sample_role=sample.sample_role,
        control_type=sample.control_type,
        strain_proxy_score=proxy_score,
        memory_proxy_recoverability=memory_score,
        event_localization_score=event_score,
        chirp_order_score=chirp_score,
        composition_alignment_score=composition,
        permanence_score=permanence,
        baseline_quiet_score=baseline_quiet,
        generic_band_energy_score=generic,
        event_energy_fraction=event_fraction,
        peak_time_error=peak_error,
        centroid_slope_hz=centroid_slope,
        step_profile_correlation=step_corr,
        plateau_stability=plateau,
    )


def build_samples(record: StrainRecord, config: ProxyConfig) -> list[SampleSpec]:
    rng = np.random.default_rng(config.seed + 17)
    target = normalize_signal(record.strain)
    return [
        SampleSpec("strain_proxy_target", "target", "none", target),
        SampleSpec("time_shuffle_control", "negative_control", "global_time_shuffle", time_shuffle(target, rng)),
        SampleSpec("phase_scramble_control", "negative_control", "fft_phase_scramble", phase_scramble(target, rng)),
        SampleSpec(
            "amplitude_preserving_surrogate_control",
            "negative_control",
            "iaaft_amplitude_spectrum_surrogate",
            amplitude_preserving_surrogate(target, rng),
        ),
        SampleSpec(
            "event_window_scramble_control",
            "negative_control",
            "event_window_chunk_scramble",
            event_window_scramble(record, config, target, rng),
        ),
        SampleSpec(
            "off_event_injection_control",
            "negative_control",
            "structured_event_shifted_out_of_window",
            off_event_shift(record, config, target),
        ),
        SampleSpec("noise_only_baseline_control", "negative_control", "noise_only_same_scale", noise_only(record, config, target, rng)),
    ]


def time_shuffle(strain: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    return strain[rng.permutation(strain.size)]


def phase_scramble(strain: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    spectrum = np.fft.rfft(normalize_signal(strain))
    phases = rng.uniform(0.0, 2.0 * math.pi, size=spectrum.size)
    phases[0] = 0.0
    if spectrum.size > 1:
        phases[-1] = 0.0
    scrambled = np.abs(spectrum) * np.exp(1j * phases)
    return normalize_signal(np.fft.irfft(scrambled, n=strain.size))


def amplitude_preserving_surrogate(strain: np.ndarray, rng: np.random.Generator, iterations: int = 24) -> np.ndarray:
    target = normalize_signal(strain)
    sorted_target = np.sort(target)
    target_amplitude = np.abs(np.fft.rfft(target))
    surrogate = rng.permutation(target)
    for _ in range(iterations):
        spectrum = np.fft.rfft(surrogate)
        surrogate = np.fft.irfft(target_amplitude * np.exp(1j * np.angle(spectrum)), n=target.size)
        ranks = np.argsort(np.argsort(surrogate))
        surrogate = sorted_target[ranks]
    return normalize_signal(surrogate)


def event_window_scramble(
    record: StrainRecord,
    config: ProxyConfig,
    strain: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.array(strain, copy=True)
    mask = event_masks(record.time, record.event_time, config)["event"]
    indices = np.flatnonzero(mask)
    if indices.size < 16:
        return time_shuffle(out, rng)
    chunk_size = max(4, int(round(0.026 * record.sample_rate)))
    chunks = [indices[start : start + chunk_size] for start in range(0, indices.size, chunk_size)]
    values = [out[chunk].copy() for chunk in chunks]
    order = rng.permutation(len(values))
    shuffled_pieces = []
    for idx in order:
        piece = values[int(idx)]
        shuffled_pieces.append(piece[::-1] if rng.random() < 0.5 else piece)
    rebuilt = np.concatenate(shuffled_pieces)
    out[indices] = rebuilt[: indices.size]
    return normalize_signal(out)


def off_event_shift(record: StrainRecord, config: ProxyConfig, strain: np.ndarray) -> np.ndarray:
    shift_seconds = max(0.32, config.duration * 0.24)
    shift_count = int(round(shift_seconds * record.sample_rate))
    return normalize_signal(np.roll(strain, shift_count))


def noise_only(record: StrainRecord, config: ProxyConfig, strain: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    masks = event_masks(record.time, record.event_time, config)
    if np.any(masks["pre"]):
        std = float(np.std(strain[masks["pre"]]))
    else:
        std = float(np.std(strain))
    return normalize_signal(rng.normal(scale=max(std, 1.0e-6), size=strain.size))


def run_probe(record: StrainRecord, config: ProxyConfig) -> tuple[list[SampleSpec], list[SampleSummary], dict[str, Any]]:
    samples = build_samples(record, config)
    summaries = [evaluate_sample(record, config, sample) for sample in samples]
    status = build_status(record, config, summaries)
    return samples, summaries, status


def build_status(record: StrainRecord, config: ProxyConfig, summaries: list[SampleSummary]) -> dict[str, Any]:
    target = next(row for row in summaries if row.sample_role == "target")
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    control_scores = np.array([row.strain_proxy_score for row in controls], dtype=float)
    control_threshold = target.strain_proxy_score - config.pass_min_margin
    false_positive_rate = float(np.mean(control_scores >= control_threshold))
    control_max = max(controls, key=lambda row: row.strain_proxy_score)
    margin = float(target.strain_proxy_score - control_max.strain_proxy_score)
    separation_z = float(
        (target.strain_proxy_score - float(np.mean(control_scores)))
        / max(float(np.std(control_scores)), 1.0e-9)
    )
    target_passes = bool(
        target.strain_proxy_score >= config.pass_min_proxy_score
        and target.memory_proxy_recoverability >= config.pass_min_memory_recoverability
        and target.event_localization_score >= config.pass_min_event_localization
        and target.chirp_order_score >= config.pass_min_chirp_order
        and target.composition_alignment_score >= config.pass_min_composition_alignment
    )
    controls_pass = bool(
        false_positive_rate <= config.pass_max_control_false_positive_rate
        and margin >= config.pass_min_margin
        and separation_z >= config.pass_min_separation_z
    )
    if target_passes and controls_pass:
        bridge_status = "PASS"
        failure_mode = "NONE"
    elif target.strain_proxy_score >= 0.62 and separation_z >= 0.85:
        bridge_status = "MARGINAL"
        failure_mode = "GW_PROXY_CONTROL_SEPARATION_INCOMPLETE"
    else:
        bridge_status = "OPEN"
        failure_mode = "STRAIN_DERIVED_PROXY_NOT_DISTINCT"

    return {
        "phase": "64",
        "bridge_name": "gw_memory_entry_gate",
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "success_definition": (
            "PASS means HAOS-style strain-derived memory/composition proxy telemetry "
            "distinguishes a structured strain-like event from time-shuffle, phase-scramble, "
            "amplitude-preserving, event-window, off-event, and noise controls. It does not "
            "mean real gravitational memory, BMS charge, soft theorem, celestial holography, "
            "or S-matrix recovery."
        ),
        "data_source": {
            "source_label": record.source_label,
            "source_kind": record.source_kind,
            "detector": record.detector,
            "sample_rate": float(record.sample_rate),
            "sample_count": int(record.strain.size),
            "event_time": float(record.event_time),
            "claim_label": "strain-derived proxy input",
        },
        "config": asdict(config),
        "target_strain_proxy_score": float(target.strain_proxy_score),
        "target_memory_proxy_recoverability": float(target.memory_proxy_recoverability),
        "target_event_localization_score": float(target.event_localization_score),
        "target_chirp_order_score": float(target.chirp_order_score),
        "target_composition_alignment_score": float(target.composition_alignment_score),
        "target_permanence_score": float(target.permanence_score),
        "control_max_score": float(control_max.strain_proxy_score),
        "control_max_sample": control_max.sample,
        "control_false_positive_threshold": float(control_threshold),
        "control_false_positive_rate": false_positive_rate,
        "score_margin_over_best_control": margin,
        "separation_z": separation_z,
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "65_real_data_replicates_or_phase64_claim_gate_review",
        "non_claims": [
            "no real gravitational-memory observable claim",
            "no BMS charge or soft-hair recovery claim",
            "no celestial holography recovery claim",
            "no real soft theorem claim",
            "no S-matrix reconstruction claim",
            "strain-derived proxy telemetry only",
        ],
    }


def write_csv(path: Path, summaries: list[SampleSummary]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(row) for row in summaries)


def write_status(path: Path, status: dict[str, Any]) -> None:
    path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")


def write_report(path: Path, record: StrainRecord, summaries: list[SampleSummary], status: dict[str, Any]) -> None:
    rows = "\n".join(
        "| {sample} | {sample_role} | {control_type} | {strain_proxy_score:.6f} | "
        "{memory_proxy_recoverability:.6f} | {event_localization_score:.6f} | "
        "{chirp_order_score:.6f} | {composition_alignment_score:.6f} |".format(**asdict(row))
        for row in summaries
    )
    target = next(row for row in summaries if row.sample_role == "target")
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    best_control = max(controls, key=lambda row: row.strain_proxy_score)
    report = f"""# Phase 64 GW Memory Entry Gate

## Scope

This report is generated by `experiments/physics_bridge/gw_memory_entry_gate/run_gw_memory_proxy.py`.

Phase 64 is a real-data entry gate for HAOS-IIP celestial-facing telemetry. It tests a strain-derived proxy on a local strain file or, by default, on a deterministic GW150914-like surrogate. The default surrogate keeps the phase reproducible without external downloads.

This is not a gravitational-memory detection. It is not a BMS charge, soft theorem, celestial holography, or S-matrix result.

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- source: `{record.source_label}`
- source kind: `{record.source_kind}`
- detector: `{record.detector}`
- target strain proxy score: `{target.strain_proxy_score:.6f}`
- best control: `{best_control.sample}` at `{best_control.strain_proxy_score:.6f}`
- margin over best control: `{status['score_margin_over_best_control']:.6f}`
- separation z: `{status['separation_z']:.6f}`
- false-positive rate: `{status['control_false_positive_rate']:.6f}`

## Diagnostics

| Sample | Role | Control | Proxy score | Memory proxy | Event localization | Chirp order | Composition alignment |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
{rows}

## Interpretation

The target row is allowed to mean only this: HAOS-style telemetry separated a structured strain-like event from the specified controls under a strain-derived proxy. If the input is the default surrogate, the result is a surrogate-data entry-gate result. If the input is a local GWOSC-style file, the result is still a strain-derived proxy result unless a separate, external gravitational-memory observable is supplied and validated.

## Controls

- `global_time_shuffle`: preserves marginal strain values while breaking temporal order.
- `fft_phase_scramble`: preserves the Fourier amplitude spectrum while breaking event phase.
- `iaaft_amplitude_spectrum_surrogate`: approximately preserves amplitude distribution and spectrum.
- `event_window_chunk_scramble`: preserves much of the event-window energy while breaking chirp order.
- `structured_event_shifted_out_of_window`: keeps the structured event but moves it away from the claimed event time.
- `noise_only_same_scale`: matched-scale noise baseline.

## Claim Gate

Allowed language:

- strain-derived memory proxy
- GW-like or local strain entry gate
- structured event separated from controls

Disallowed language:

- real gravitational memory detected
- BMS charge or soft hair recovered
- celestial holography recovered
- soft theorem recovered
- S-matrix reconstructed

## Generated Artifacts

- CSV: `experiments/physics_bridge/gw_memory_entry_gate/gw_memory_diagnostics.csv`
- JSON: `experiments/physics_bridge/gw_memory_entry_gate/bridge_status.json`
- prepared strain preview: `experiments/physics_bridge/gw_memory_entry_gate/prepared_strain_preview.csv`
- prepared strain metadata: `experiments/physics_bridge/gw_memory_entry_gate/prepared_strain_metadata.json`
- figures: `experiments/physics_bridge/gw_memory_entry_gate/figures/`
"""
    path.write_text(report, encoding="utf-8")


def write_figures(output_dir: Path, record: StrainRecord, samples: list[SampleSpec], summaries: list[SampleSummary], config: ProxyConfig) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        (output_dir / "figures_not_generated.txt").write_text(
            "Matplotlib is not installed; rerun with `uv run --with numpy --with matplotlib ...`.\n",
            encoding="utf-8",
        )
        return

    figure_dir = output_dir / FIGURE_DIR.name
    figure_dir.mkdir(parents=True, exist_ok=True)
    target_sample = next(sample for sample in samples if sample.sample_role == "target")
    target_products = compute_products(record, config, target_sample.strain)

    plt.figure(figsize=(10.0, 7.0))
    ax1 = plt.subplot(3, 1, 1)
    ax1.plot(record.time, normalize_signal(target_sample.strain), lw=0.8, color="#273043")
    ax1.axvline(record.event_time, color="#bf4342", lw=1.2)
    ax1.set_ylabel("strain z")
    ax1.set_title("Prepared strain-derived proxy input")
    ax2 = plt.subplot(3, 1, 2, sharex=ax1)
    ax2.plot(record.time, target_products.envelope, lw=1.0, color="#2f7f95", label="envelope")
    ax2.plot(record.time, target_products.flux / max(float(np.max(target_products.flux)), 1.0e-12), lw=1.0, color="#e0a458", label="normalized flux")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.set_ylabel("event energy")
    ax3 = plt.subplot(3, 1, 3, sharex=ax1)
    ax3.plot(record.time, target_products.cumulative_memory_proxy, lw=1.2, color="#335c67", label="memory proxy")
    ax3.plot(record.time, target_products.ideal_gate, lw=1.0, color="#9e2a2b", ls="--", label="event gate")
    ax3.legend(loc="lower right", fontsize=8)
    ax3.set_xlabel("time (s)")
    ax3.set_ylabel("cumulative proxy")
    plt.tight_layout()
    plt.savefig(figure_dir / "strain_proxy_trace.png", dpi=180)
    plt.close()

    ordered = sorted(summaries, key=lambda row: row.strain_proxy_score, reverse=True)
    colors = ["#1f7a5c" if row.sample_role == "target" else "#8a8f98" for row in ordered]
    plt.figure(figsize=(10.0, 4.8))
    plt.bar([row.sample for row in ordered], [row.strain_proxy_score for row in ordered], color=colors)
    plt.xticks(rotation=32, ha="right", fontsize=8)
    plt.ylabel("strain proxy score")
    plt.title("Phase 64 target/control separation")
    plt.tight_layout()
    plt.savefig(figure_dir / "control_comparison.png", dpi=180)
    plt.close()

    target = next(row for row in summaries if row.sample_role == "target")
    best_control = max((row for row in summaries if row.sample_role == "negative_control"), key=lambda row: row.strain_proxy_score)
    labels = ["memory", "event", "chirp", "composition", "permanence"]
    target_values = [
        target.memory_proxy_recoverability,
        target.event_localization_score,
        target.chirp_order_score,
        target.composition_alignment_score,
        target.permanence_score,
    ]
    control_values = [
        best_control.memory_proxy_recoverability,
        best_control.event_localization_score,
        best_control.chirp_order_score,
        best_control.composition_alignment_score,
        best_control.permanence_score,
    ]
    x = np.arange(len(labels))
    plt.figure(figsize=(8.2, 4.8))
    plt.bar(x - 0.18, target_values, width=0.36, label="target", color="#1f7a5c")
    plt.bar(x + 0.18, control_values, width=0.36, label=best_control.sample, color="#8a8f98")
    plt.xticks(x, labels)
    plt.ylim(0.0, 1.05)
    plt.ylabel("score")
    plt.title("Metric breakdown")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figure_dir / "metric_breakdown.png", dpi=180)
    plt.close()

    spectrogram, freqs, centers = short_time_spectrum(target_products.filtered, record.sample_rate, config)
    plt.figure(figsize=(9.2, 4.8))
    plt.imshow(
        spectrogram,
        aspect="auto",
        origin="lower",
        extent=[centers[0], centers[-1], freqs[0], freqs[-1]],
        cmap="magma",
    )
    plt.axvline(record.event_time, color="white", lw=1.0)
    plt.colorbar(label="normalized power")
    plt.xlabel("time (s)")
    plt.ylabel("frequency (Hz)")
    plt.title("Strain-derived proxy time-frequency view")
    plt.tight_layout()
    plt.savefig(figure_dir / "time_frequency_proxy.png", dpi=180)
    plt.close()


def short_time_spectrum(
    filtered: np.ndarray,
    sample_rate: float,
    config: ProxyConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    window_count = max(32, int(round(0.12 * sample_rate)))
    hop = max(8, window_count // 4)
    spectra: list[np.ndarray] = []
    centers: list[float] = []
    for start in range(0, max(1, filtered.size - window_count), hop):
        segment = filtered[start : start + window_count]
        if segment.size < window_count:
            break
        window = np.hanning(segment.size)
        power = np.abs(np.fft.rfft(segment * window)) ** 2
        freqs = np.fft.rfftfreq(segment.size, d=1.0 / float(sample_rate))
        mask = (freqs >= config.f_min) & (freqs <= config.f_max)
        spectra.append(power[mask])
        centers.append((start + window_count / 2.0) / float(sample_rate))
    if not spectra:
        return np.zeros((1, 1)), np.array([0.0]), np.array([0.0])
    data = np.array(spectra, dtype=float).T
    data = data / max(float(np.max(data)), 1.0e-12)
    freqs = np.fft.rfftfreq(window_count, d=1.0 / float(sample_rate))
    mask = (freqs >= config.f_min) & (freqs <= config.f_max)
    return data, freqs[mask], np.array(centers, dtype=float)


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(resolved)


def main() -> None:
    args = parse_args()
    config = make_config(args)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / FIGURE_DIR.name).mkdir(parents=True, exist_ok=True)

    raw_record = load_strain_record(
        args.input_file,
        sample_rate=float(args.sample_rate),
        duration=config.duration,
        event_time=config.event_time,
        seed=config.seed,
        signal_to_noise=config.signal_to_noise,
        detector=str(args.detector),
    )
    record = prepare_record(
        raw_record,
        duration=config.duration,
        analysis_sample_rate=config.analysis_sample_rate,
        event_time=config.event_time,
    )

    samples, summaries, status = run_probe(record, config)
    write_csv(output_dir / CSV_PATH.name, summaries)
    write_status(output_dir / STATUS_PATH.name, status)
    write_report(output_dir / REPORT_PATH.name, record, summaries, status)
    write_record_preview(output_dir / PREVIEW_PATH.name, record)
    write_record_metadata(output_dir / METADATA_PATH.name, record)
    write_figures(output_dir, record, samples, summaries, config)

    print("Phase 64 GW memory entry gate generated.")
    print(f"CSV: {relative(output_dir / CSV_PATH.name)}")
    print(f"JSON: {relative(output_dir / STATUS_PATH.name)}")
    print(f"Report: {relative(output_dir / REPORT_PATH.name)}")
    print(f"bridge_status: {status['bridge_status']} ({status['failure_mode']})")
    print(
        f"target_score={status['target_strain_proxy_score']:.6f}; "
        f"control_max={status['control_max_score']:.6f}; "
        f"margin={status['score_margin_over_best_control']:.6f}; "
        f"z={status['separation_z']:.6f}"
    )


if __name__ == "__main__":
    main()
