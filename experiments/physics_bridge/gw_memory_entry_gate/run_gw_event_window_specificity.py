#!/usr/bin/env python3
"""Phase 66 event-internal specificity for the GW entry gate.

Phase 65 found that envelope-locked and sliding/micro event-window controls can
still compete with the strain-derived proxy. Phase 66 adds a narrower
event-internal diagnostic: does the candidate event preserve a coherent
chirp-like time-frequency ridge and phase-continuity structure, rather than only
the event timing envelope?

This is still a claim-gated proxy. A PASS here means event-window specificity
improved on deterministic strain-like surrogate replicates. It does not mean
real gravitational memory, BMS charges, soft theorems, celestial holography, or
S-matrix recovery.
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

from run_gw_event_window_hardening import (
    HardeningConfig,
    build_hardening_samples,
    build_records,
    event_indices,
    proxy_config,
)
from run_gw_memory_proxy import (
    ProxyConfig,
    SampleSummary,
    evaluate_sample,
    fft_bandpass,
    moving_average,
    normalize_signal,
    pearson_correlation,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(__file__).resolve().parent

DIAGNOSTICS_PATH = ROOT / "event_window_specificity_diagnostics.csv"
SUMMARY_PATH = ROOT / "event_window_specificity_summary.csv"
STATUS_PATH = ROOT / "event_window_specificity_status.json"
REPORT_PATH = ROOT / "gw_event_window_specificity_report.md"
FIGURE_DIR = ROOT / "figures"

DIAGNOSTIC_FIELDNAMES = [
    "event_id",
    "seed",
    "signal_to_noise",
    "sample",
    "sample_role",
    "control_type",
    "specificity_score",
    "base_strain_proxy_score",
    "ridge_coherence_score",
    "phase_continuity_score",
    "ridge_monotonicity",
    "ridge_slope_hz",
    "ridge_concentration",
    "ridge_smoothness",
    "phase_valid_fraction",
    "phase_frequency_correlation",
    "phase_smoothness",
    "phase_range_hz",
]

SUMMARY_FIELDNAMES = [
    "event_id",
    "target_specificity_score",
    "best_control_sample",
    "best_control_type",
    "best_control_specificity_score",
    "specificity_margin",
    "event_window_specificity_margin",
    "separation_z",
    "control_false_positive_rate",
]


@dataclass(frozen=True)
class SpecificityConfig:
    event_count: int = 12
    start_seed: int = 6601
    duration: float = 2.0
    sample_rate: float = 1024.0
    analysis_sample_rate: float = 1024.0
    signal_to_noise: float = 14.0
    signal_to_noise_jitter: float = 3.5
    pass_min_mean_target_score: float = 0.78
    pass_min_mean_margin: float = 0.25
    pass_min_min_margin: float = 0.12
    pass_min_mean_event_window_margin: float = 0.24
    pass_min_mean_separation_z: float = 1.50
    pass_max_control_false_positive_rate: float = 0.0


@dataclass(frozen=True)
class RidgeDiagnostics:
    score: float
    monotonicity: float
    slope_hz: float
    concentration: float
    smoothness: float


@dataclass(frozen=True)
class PhaseDiagnostics:
    score: float
    valid_fraction: float
    frequency_correlation: float
    smoothness: float
    range_hz: float


@dataclass(frozen=True)
class SpecificitySummary:
    event_id: str
    seed: int
    signal_to_noise: float
    sample: str
    sample_role: str
    control_type: str
    specificity_score: float
    base_strain_proxy_score: float
    ridge_coherence_score: float
    phase_continuity_score: float
    ridge_monotonicity: float
    ridge_slope_hz: float
    ridge_concentration: float
    ridge_smoothness: float
    phase_valid_fraction: float
    phase_frequency_correlation: float
    phase_smoothness: float
    phase_range_hz: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 66 GW event-window specificity.")
    parser.add_argument("--event-count", type=int, default=12)
    parser.add_argument("--start-seed", type=int, default=6601)
    parser.add_argument("--duration", type=float, default=2.0)
    parser.add_argument("--sample-rate", type=float, default=1024.0)
    parser.add_argument("--analysis-sample-rate", type=float, default=1024.0)
    parser.add_argument("--signal-to-noise", type=float, default=14.0)
    parser.add_argument("--signal-to-noise-jitter", type=float, default=3.5)
    parser.add_argument("--input-files", type=Path, nargs="*", default=[])
    parser.add_argument("--event-times", type=str, default="")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated CSV, JSON, report, and figures.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> SpecificityConfig:
    return SpecificityConfig(
        event_count=int(args.event_count),
        start_seed=int(args.start_seed),
        duration=float(args.duration),
        sample_rate=float(args.sample_rate),
        analysis_sample_rate=float(args.analysis_sample_rate),
        signal_to_noise=float(args.signal_to_noise),
        signal_to_noise_jitter=float(args.signal_to_noise_jitter),
    )


def hardening_config(config: SpecificityConfig) -> HardeningConfig:
    return HardeningConfig(
        event_count=config.event_count,
        start_seed=config.start_seed,
        duration=config.duration,
        sample_rate=config.sample_rate,
        analysis_sample_rate=config.analysis_sample_rate,
        signal_to_noise=config.signal_to_noise,
        signal_to_noise_jitter=config.signal_to_noise_jitter,
    )


def analytic_signal(values: np.ndarray) -> np.ndarray:
    """Small NumPy-only Hilbert transform helper."""

    n = values.size
    spectrum = np.fft.fft(values)
    weights = np.zeros(n, dtype=float)
    if n % 2 == 0:
        weights[0] = 1.0
        weights[n // 2] = 1.0
        weights[1 : n // 2] = 2.0
    else:
        weights[0] = 1.0
        weights[1 : (n + 1) // 2] = 2.0
    return np.fft.ifft(spectrum * weights)


def event_windows(record_time: np.ndarray, center: float, indices: np.ndarray, sample_rate: float) -> list[np.ndarray]:
    if indices.size < 32:
        return []
    window = max(24, int(round(0.052 * sample_rate)))
    hop = max(6, int(round(0.014 * sample_rate)))
    windows: list[np.ndarray] = []
    local = indices
    for start in range(0, max(1, local.size - window), hop):
        item = local[start : start + window]
        if item.size == window:
            windows.append(item)
    if not windows:
        windows.append(local)
    return windows


def ridge_coherence_for_record(record: Any, config: ProxyConfig, strain: np.ndarray) -> RidgeDiagnostics:
    filtered = fft_bandpass(strain, record.sample_rate, config.f_min, config.f_max)
    indices = event_indices(record, config)
    windows = event_windows(record.time, record.event_time, indices, record.sample_rate)
    if len(windows) < 3:
        return RidgeDiagnostics(0.0, 0.0, 0.0, 0.0, 0.0)

    centroids: list[float] = []
    concentrations: list[float] = []
    for item in windows:
        segment = normalize_signal(filtered[item])
        window = np.hanning(segment.size)
        power = np.abs(np.fft.rfft(segment * window)) ** 2
        freqs = np.fft.rfftfreq(segment.size, d=1.0 / float(record.sample_rate))
        mask = (freqs >= config.f_min) & (freqs <= config.f_max)
        power = power[mask]
        freqs = freqs[mask]
        total = float(np.sum(power))
        if total <= 1.0e-12:
            continue
        centroid = float(np.sum(freqs * power) / total)
        peak_idx = int(np.argmax(power))
        band = np.abs(freqs - freqs[peak_idx]) <= 28.0
        concentration = float(np.sum(power[band]) / total)
        centroids.append(centroid)
        concentrations.append(concentration)

    if len(centroids) < 3:
        return RidgeDiagnostics(0.0, 0.0, 0.0, 0.0, 0.0)

    center_values = np.asarray(centroids, dtype=float)
    diffs = np.diff(center_values)
    positive = diffs > 0.8
    monotonicity = float(np.mean(positive))
    slope_hz = float(center_values[-1] - center_values[0])
    slope_score = float(np.clip(slope_hz / 82.0, 0.0, 1.0))
    concentration = float(np.mean(concentrations))
    concentration_score = float(np.clip((concentration - 0.30) / 0.46, 0.0, 1.0))
    if diffs.size >= 2:
        curvature = float(np.median(np.abs(np.diff(diffs))))
        step = float(np.median(np.abs(diffs))) + 1.0
        smoothness = float(np.clip(1.0 - curvature / (step + 16.0), 0.0, 1.0))
    else:
        smoothness = 0.0
    score = float(0.34 * monotonicity + 0.28 * slope_score + 0.24 * concentration_score + 0.14 * smoothness)
    return RidgeDiagnostics(score, monotonicity, slope_hz, concentration, smoothness)


def phase_continuity_for_record(record: Any, config: ProxyConfig, strain: np.ndarray) -> PhaseDiagnostics:
    filtered = fft_bandpass(strain, record.sample_rate, config.f_min, config.f_max)
    analytic = analytic_signal(filtered)
    phase = np.unwrap(np.angle(analytic))
    dt = 1.0 / float(record.sample_rate)
    inst_freq = np.gradient(phase, dt) / (2.0 * math.pi)
    smooth_count = max(3, int(round(0.018 * record.sample_rate)))
    inst_freq = moving_average(inst_freq, smooth_count)
    indices = event_indices(record, config)
    if indices.size < 24:
        return PhaseDiagnostics(0.0, 0.0, 0.0, 0.0, 0.0)

    freq = inst_freq[indices]
    time = record.time[indices]
    finite = np.isfinite(freq)
    valid_band = finite & (freq >= config.f_min * 0.60) & (freq <= config.f_max * 1.35)
    valid_fraction = float(np.mean(valid_band))
    if np.count_nonzero(valid_band) < 8:
        return PhaseDiagnostics(0.0, valid_fraction, 0.0, 0.0, 0.0)

    freq = freq[valid_band]
    time = time[valid_band]
    corr = float(np.clip(max(pearson_correlation(freq, time), 0.0), 0.0, 1.0))
    freq_range = float(np.percentile(freq, 88) - np.percentile(freq, 12))
    range_score = float(np.clip(freq_range / 155.0, 0.0, 1.0))
    roughness = float(np.median(np.abs(np.diff(freq)))) if freq.size > 2 else 1.0e9
    smoothness = float(np.clip(1.0 - roughness / 72.0, 0.0, 1.0))
    score = float(0.36 * corr + 0.24 * smoothness + 0.22 * valid_fraction + 0.18 * range_score)
    return PhaseDiagnostics(score, valid_fraction, corr, smoothness, freq_range)


def evaluate_specificity(record: Any, config: ProxyConfig, seed: int, event_id: str) -> list[SpecificitySummary]:
    samples = build_hardening_samples(record, config, seed)
    rows: list[SpecificitySummary] = []
    for sample in samples:
        base = evaluate_sample(record, config, sample)
        ridge = ridge_coherence_for_record(record, config, sample.strain)
        phase = phase_continuity_for_record(record, config, sample.strain)
        specificity = float(
            0.34 * base.strain_proxy_score
            + 0.30 * ridge.score
            + 0.27 * phase.score
            + 0.09 * base.composition_alignment_score
        )
        rows.append(
            SpecificitySummary(
                event_id=event_id,
                seed=seed,
                signal_to_noise=float(config.signal_to_noise),
                sample=sample.sample,
                sample_role=sample.sample_role,
                control_type=sample.control_type,
                specificity_score=specificity,
                base_strain_proxy_score=base.strain_proxy_score,
                ridge_coherence_score=ridge.score,
                phase_continuity_score=phase.score,
                ridge_monotonicity=ridge.monotonicity,
                ridge_slope_hz=ridge.slope_hz,
                ridge_concentration=ridge.concentration,
                ridge_smoothness=ridge.smoothness,
                phase_valid_fraction=phase.valid_fraction,
                phase_frequency_correlation=phase.frequency_correlation,
                phase_smoothness=phase.smoothness,
                phase_range_hz=phase.range_hz,
            )
        )
    return rows


def summarize_event(rows: list[SpecificitySummary], config: ProxyConfig) -> dict[str, Any]:
    target = next(row for row in rows if row.sample_role == "target")
    controls = [row for row in rows if row.sample_role == "negative_control"]
    control_scores = np.array([row.specificity_score for row in controls], dtype=float)
    threshold = target.specificity_score - config.pass_min_margin
    best = max(controls, key=lambda row: row.specificity_score)
    event_controls = [
        row
        for row in controls
        if "event" in row.control_type or "chirp" in row.control_type or "timing" in row.control_type
    ]
    event_best = max(event_controls, key=lambda row: row.specificity_score)
    return {
        "event_id": target.event_id,
        "target_specificity_score": float(target.specificity_score),
        "best_control_sample": best.sample,
        "best_control_type": best.control_type,
        "best_control_specificity_score": float(best.specificity_score),
        "specificity_margin": float(target.specificity_score - best.specificity_score),
        "event_window_specificity_margin": float(target.specificity_score - event_best.specificity_score),
        "separation_z": float(
            (target.specificity_score - float(np.mean(control_scores)))
            / max(float(np.std(control_scores)), 1.0e-9)
        ),
        "control_false_positive_rate": float(np.mean(control_scores >= threshold)),
    }


def run_specificity(args: argparse.Namespace, config: SpecificityConfig) -> tuple[list[SpecificitySummary], list[dict[str, Any]], dict[str, Any]]:
    hconfig = hardening_config(config)
    all_rows: list[SpecificitySummary] = []
    event_rows: list[dict[str, Any]] = []
    for event_id, seed, snr, record in build_records(args, hconfig):
        pconfig = proxy_config(hconfig, seed, snr)
        rows = evaluate_specificity(record, pconfig, seed, event_id)
        all_rows.extend(rows)
        event_rows.append(summarize_event(rows, pconfig))
    status = build_status(config, event_rows)
    return all_rows, event_rows, status


def build_status(config: SpecificityConfig, event_rows: list[dict[str, Any]]) -> dict[str, Any]:
    margins = np.array([row["specificity_margin"] for row in event_rows], dtype=float)
    event_margins = np.array([row["event_window_specificity_margin"] for row in event_rows], dtype=float)
    target_scores = np.array([row["target_specificity_score"] for row in event_rows], dtype=float)
    z_values = np.array([row["separation_z"] for row in event_rows], dtype=float)
    fprs = np.array([row["control_false_positive_rate"] for row in event_rows], dtype=float)
    weakest = min(event_rows, key=lambda row: row["specificity_margin"])
    highest = max(event_rows, key=lambda row: row["best_control_specificity_score"])
    aggregate_fpr = float(np.mean(fprs > 0.0))

    target_passes = bool(float(np.mean(target_scores)) >= config.pass_min_mean_target_score)
    controls_pass = bool(
        float(np.mean(margins)) >= config.pass_min_mean_margin
        and float(np.min(margins)) >= config.pass_min_min_margin
        and float(np.mean(event_margins)) >= config.pass_min_mean_event_window_margin
        and float(np.mean(z_values)) >= config.pass_min_mean_separation_z
        and aggregate_fpr <= config.pass_max_control_false_positive_rate
    )
    if target_passes and controls_pass:
        bridge_status = "PASS"
        failure_mode = "NONE"
    elif target_passes and float(np.mean(z_values)) >= 1.0:
        bridge_status = "MARGINAL"
        failure_mode = "EVENT_INTERNAL_SPECIFICITY_PARTIAL"
    else:
        bridge_status = "OPEN"
        failure_mode = "EVENT_INTERNAL_SPECIFICITY_NOT_DISTINCT"

    return {
        "phase": "66",
        "bridge_name": "gw_event_window_specificity",
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "success_definition": (
            "PASS means event-internal ridge coherence and phase-continuity metrics "
            "separate structured strain-like surrogate events from envelope-preserving "
            "event-window controls across replicates. It does not mean real gravitational "
            "memory, BMS charges, soft theorems, celestial holography, or S-matrix recovery."
        ),
        "config": asdict(config),
        "event_count": len(event_rows),
        "mean_target_specificity_score": float(np.mean(target_scores)),
        "std_target_specificity_score": float(np.std(target_scores)),
        "mean_specificity_margin": float(np.mean(margins)),
        "min_specificity_margin": float(np.min(margins)),
        "mean_event_window_specificity_margin": float(np.mean(event_margins)),
        "min_event_window_specificity_margin": float(np.min(event_margins)),
        "mean_separation_z": float(np.mean(z_values)),
        "min_separation_z": float(np.min(z_values)),
        "aggregate_control_false_positive_rate": aggregate_fpr,
        "weakest_event_id": str(weakest["event_id"]),
        "weakest_event_margin": float(weakest["specificity_margin"]),
        "highest_control_event_id": str(highest["event_id"]),
        "highest_control_sample": str(highest["best_control_sample"]),
        "highest_control_score": float(highest["best_control_specificity_score"]),
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "67_real_gwosc_replication_or_event_window_metric_ablation",
        "non_claims": [
            "no real gravitational-memory observable claim",
            "no BMS charge or soft-hair recovery claim",
            "no celestial holography recovery claim",
            "no real soft theorem claim",
            "no S-matrix reconstruction claim",
            "event-internal strain-derived proxy telemetry only",
        ],
    }


def write_diagnostics(path: Path, rows: list[SpecificitySummary]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DIAGNOSTIC_FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in SUMMARY_FIELDNAMES})


def write_status(path: Path, status: dict[str, Any]) -> None:
    path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")


def write_report(path: Path, rows: list[dict[str, Any]], status: dict[str, Any]) -> None:
    table_rows = "\n".join(
        "| {event_id} | {target_specificity_score:.6f} | {best_control_sample} | "
        "{best_control_specificity_score:.6f} | {specificity_margin:.6f} | "
        "{event_window_specificity_margin:.6f} | {separation_z:.6f} |".format(**row)
        for row in rows
    )
    report = f"""# Phase 66 GW Event-Window Specificity

## Scope

This report is generated by `experiments/physics_bridge/gw_memory_entry_gate/run_gw_event_window_specificity.py`.

Phase 66 responds to the Phase 65 `EVENT_WINDOW_CONTROL_LEAK_REMAINS` result by adding event-internal ridge-coherence and phase-continuity metrics. The goal is to distinguish a coherent chirp-like event interior from controls that preserve only timing-envelope structure.

This is not a gravitational-memory detection. It is not a BMS charge, soft theorem, celestial holography, or S-matrix result.

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- event count: `{status['event_count']}`
- mean target specificity: `{status['mean_target_specificity_score']:.6f}`
- mean specificity margin: `{status['mean_specificity_margin']:.6f}`
- minimum specificity margin: `{status['min_specificity_margin']:.6f}`
- mean event-window specificity margin: `{status['mean_event_window_specificity_margin']:.6f}`
- mean separation z: `{status['mean_separation_z']:.6f}`
- aggregate false-positive rate: `{status['aggregate_control_false_positive_rate']:.6f}`
- weakest event: `{status['weakest_event_id']}` at margin `{status['weakest_event_margin']:.6f}`
- highest control: `{status['highest_control_sample']}` in `{status['highest_control_event_id']}` at `{status['highest_control_score']:.6f}`

## Event Summary

| Event | Target specificity | Best control | Control specificity | Margin | Event-window margin | z |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
{table_rows}

## Added Metrics

- `ridge_coherence_score`: checks event-window time-frequency ridge monotonicity, slope, concentration, and smoothness.
- `phase_continuity_score`: checks event-window analytic-phase frequency validity, temporal correlation, smoothness, and range.

## Interpretation

Allowed language:

- event-internal strain-derived specificity
- chirp/ridge/phase-continuity proxy
- envelope-preserving event-window controls separated or characterized

Disallowed language:

- real gravitational memory detected
- BMS charge or soft hair recovered
- celestial holography recovered
- soft theorem recovered
- S-matrix reconstructed

## Generated Artifacts

- diagnostics CSV: `experiments/physics_bridge/gw_memory_entry_gate/event_window_specificity_diagnostics.csv`
- summary CSV: `experiments/physics_bridge/gw_memory_entry_gate/event_window_specificity_summary.csv`
- JSON: `experiments/physics_bridge/gw_memory_entry_gate/event_window_specificity_status.json`
- figures: `experiments/physics_bridge/gw_memory_entry_gate/figures/phase66_*.png`
"""
    path.write_text(report, encoding="utf-8")


def write_figures(output_dir: Path, rows: list[SpecificitySummary], event_rows: list[dict[str, Any]]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        (output_dir / "phase66_figures_not_generated.txt").write_text(
            "Matplotlib is not installed; rerun with `uv run --with numpy --with matplotlib ...`.\n",
            encoding="utf-8",
        )
        return

    figure_dir = output_dir / FIGURE_DIR.name
    figure_dir.mkdir(parents=True, exist_ok=True)
    event_ids = [row["event_id"] for row in event_rows]
    target_scores = np.array([row["target_specificity_score"] for row in event_rows], dtype=float)
    control_scores = np.array([row["best_control_specificity_score"] for row in event_rows], dtype=float)
    margins = np.array([row["specificity_margin"] for row in event_rows], dtype=float)
    x = np.arange(len(event_rows))

    plt.figure(figsize=(10.5, 4.8))
    plt.plot(x, target_scores, marker="o", color="#1f7a5c", label="target specificity")
    plt.plot(x, control_scores, marker="o", color="#a23e48", label="best control")
    plt.xticks(x, event_ids, rotation=35, ha="right", fontsize=8)
    plt.ylabel("specificity score")
    plt.title("Phase 66 event-internal specificity")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figure_dir / "phase66_specificity_scores.png", dpi=180)
    plt.close()

    plt.figure(figsize=(9.0, 4.6))
    plt.bar(event_ids, margins, color="#335c67")
    plt.axhline(0.25, color="#9e2a2b", ls="--", lw=1.0, label="mean PASS bar")
    plt.xticks(rotation=35, ha="right", fontsize=8)
    plt.ylabel("target-control margin")
    plt.title("Phase 66 specificity margins")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figure_dir / "phase66_specificity_margins.png", dpi=180)
    plt.close()

    targets = [row for row in rows if row.sample_role == "target"]
    controls = [row for row in rows if row.sample_role == "negative_control"]
    plt.figure(figsize=(7.0, 5.2))
    plt.scatter(
        [row.ridge_coherence_score for row in controls],
        [row.phase_continuity_score for row in controls],
        color="#8a8f98",
        alpha=0.75,
        label="controls",
    )
    plt.scatter(
        [row.ridge_coherence_score for row in targets],
        [row.phase_continuity_score for row in targets],
        color="#1f7a5c",
        edgecolor="black",
        linewidth=0.5,
        label="targets",
    )
    plt.xlabel("ridge coherence")
    plt.ylabel("phase continuity")
    plt.title("Phase 66 event-internal metric plane")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figure_dir / "phase66_metric_plane.png", dpi=180)
    plt.close()


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

    rows, event_rows, status = run_specificity(args, config)
    write_diagnostics(output_dir / DIAGNOSTICS_PATH.name, rows)
    write_summary(output_dir / SUMMARY_PATH.name, event_rows)
    write_status(output_dir / STATUS_PATH.name, status)
    write_report(output_dir / REPORT_PATH.name, event_rows, status)
    write_figures(output_dir, rows, event_rows)

    print("Phase 66 GW event-window specificity generated.")
    print(f"Diagnostics: {relative(output_dir / DIAGNOSTICS_PATH.name)}")
    print(f"Summary: {relative(output_dir / SUMMARY_PATH.name)}")
    print(f"JSON: {relative(output_dir / STATUS_PATH.name)}")
    print(f"Report: {relative(output_dir / REPORT_PATH.name)}")
    print(f"bridge_status: {status['bridge_status']} ({status['failure_mode']})")
    print(
        f"mean_target={status['mean_target_specificity_score']:.6f}; "
        f"mean_margin={status['mean_specificity_margin']:.6f}; "
        f"min_margin={status['min_specificity_margin']:.6f}; "
        f"mean_z={status['mean_separation_z']:.6f}"
    )


if __name__ == "__main__":
    main()
