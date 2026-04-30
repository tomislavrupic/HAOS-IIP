#!/usr/bin/env python3
"""Phase 65 event-window and multi-event hardening for the GW entry gate.

Phase 64 identified the strongest leakage path: event-window controls can
preserve timing-envelope structure while scrambling the event interior. This
runner attacks that specific leak with a stricter control ladder and multiple
deterministic strain-like surrogate events.

This remains claim-gated. A PASS here means the strain-derived proxy survives a
harder event-window control ladder across replicates. It is not a real
gravitational-memory, BMS, soft-theorem, celestial-holography, or S-matrix
claim.
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

from gw_data_loader import StrainRecord, load_strain_record, prepare_record
from run_gw_memory_proxy import (
    ProxyConfig,
    SampleSpec,
    SampleSummary,
    amplitude_preserving_surrogate,
    compute_products,
    evaluate_sample,
    event_masks,
    fft_bandpass,
    normalize_signal,
    phase_scramble,
    time_shuffle,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(__file__).resolve().parent

DIAGNOSTICS_PATH = ROOT / "event_window_hardening_diagnostics.csv"
SUMMARY_PATH = ROOT / "event_window_hardening_summary.csv"
STATUS_PATH = ROOT / "event_window_hardening_status.json"
REPORT_PATH = ROOT / "gw_event_window_hardening_report.md"
FIGURE_DIR = ROOT / "figures"

DIAGNOSTIC_FIELDNAMES = [
    "event_id",
    "source_label",
    "source_kind",
    "seed",
    "signal_to_noise",
    "sample",
    "sample_role",
    "control_type",
    "strain_proxy_score",
    "memory_proxy_recoverability",
    "event_localization_score",
    "chirp_order_score",
    "composition_alignment_score",
    "permanence_score",
    "event_energy_fraction",
    "peak_time_error",
    "centroid_slope_hz",
]

SUMMARY_FIELDNAMES = [
    "event_id",
    "source_label",
    "target_score",
    "best_control_sample",
    "best_control_type",
    "best_control_score",
    "margin_over_best_control",
    "event_window_control_max",
    "event_window_margin",
    "separation_z",
    "control_false_positive_rate",
]


@dataclass(frozen=True)
class HardeningConfig:
    event_count: int = 12
    start_seed: int = 6501
    duration: float = 2.0
    sample_rate: float = 1024.0
    analysis_sample_rate: float = 1024.0
    signal_to_noise: float = 14.0
    signal_to_noise_jitter: float = 3.5
    pass_min_mean_target_score: float = 0.78
    pass_min_mean_margin: float = 0.25
    pass_min_min_margin: float = 0.12
    pass_min_mean_event_window_margin: float = 0.20
    pass_min_mean_separation_z: float = 1.50
    pass_max_control_false_positive_rate: float = 0.0


@dataclass(frozen=True)
class EventRun:
    event_id: str
    seed: int
    signal_to_noise: float
    record: StrainRecord
    summaries: list[SampleSummary]
    summary: dict[str, Any]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 65 GW event-window hardening.")
    parser.add_argument("--event-count", type=int, default=12, help="Number of deterministic surrogate events.")
    parser.add_argument("--start-seed", type=int, default=6501)
    parser.add_argument("--duration", type=float, default=2.0)
    parser.add_argument("--sample-rate", type=float, default=1024.0)
    parser.add_argument("--analysis-sample-rate", type=float, default=1024.0)
    parser.add_argument("--signal-to-noise", type=float, default=14.0)
    parser.add_argument("--signal-to-noise-jitter", type=float, default=3.5)
    parser.add_argument(
        "--input-files",
        type=Path,
        nargs="*",
        default=[],
        help="Optional local strain files. If supplied, these replace generated surrogates.",
    )
    parser.add_argument(
        "--event-times",
        type=str,
        default="",
        help="Optional comma-separated event times matching --input-files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated CSV, JSON, report, and figures.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> HardeningConfig:
    return HardeningConfig(
        event_count=int(args.event_count),
        start_seed=int(args.start_seed),
        duration=float(args.duration),
        sample_rate=float(args.sample_rate),
        analysis_sample_rate=float(args.analysis_sample_rate),
        signal_to_noise=float(args.signal_to_noise),
        signal_to_noise_jitter=float(args.signal_to_noise_jitter),
    )


def proxy_config(config: HardeningConfig, seed: int, snr: float) -> ProxyConfig:
    return ProxyConfig(
        duration=config.duration,
        analysis_sample_rate=config.analysis_sample_rate,
        seed=seed,
        signal_to_noise=snr,
        pass_min_margin=0.14,
    )


def parse_event_times(raw: str) -> list[float | None]:
    if not raw.strip():
        return []
    values: list[float | None] = []
    for part in raw.split(","):
        text = part.strip()
        values.append(None if not text else float(text))
    return values


def build_records(args: argparse.Namespace, config: HardeningConfig) -> list[tuple[str, int, float, StrainRecord]]:
    event_times = parse_event_times(str(args.event_times))
    records: list[tuple[str, int, float, StrainRecord]] = []
    if args.input_files:
        for idx, path in enumerate(args.input_files):
            seed = config.start_seed + idx
            event_time = event_times[idx] if idx < len(event_times) else None
            raw = load_strain_record(
                path,
                sample_rate=config.sample_rate,
                duration=config.duration,
                event_time=event_time,
                seed=seed,
                signal_to_noise=config.signal_to_noise,
                detector=f"LOCAL{idx + 1}",
            )
            prepared = prepare_record(
                raw,
                duration=config.duration,
                analysis_sample_rate=config.analysis_sample_rate,
                event_time=event_time,
            )
            records.append((f"local_{idx + 1:02d}", seed, config.signal_to_noise, prepared))
        return records

    for idx in range(config.event_count):
        seed = config.start_seed + idx
        snr = float(
            np.clip(
                config.signal_to_noise + config.signal_to_noise_jitter * math.sin(1.37 * idx + 0.20),
                7.5,
                24.0,
            )
        )
        raw = load_strain_record(
            None,
            sample_rate=config.sample_rate,
            duration=config.duration,
            event_time=None,
            seed=seed,
            signal_to_noise=snr,
            detector=f"SYNTH{idx + 1:02d}",
        )
        prepared = prepare_record(
            raw,
            duration=config.duration,
            analysis_sample_rate=config.analysis_sample_rate,
            event_time=None,
        )
        records.append((f"surrogate_{idx + 1:02d}", seed, snr, prepared))
    return records


def build_hardening_samples(
    record: StrainRecord,
    config: ProxyConfig,
    seed: int,
) -> list[SampleSpec]:
    rng = np.random.default_rng(seed + 510)
    target = normalize_signal(record.strain)
    return [
        SampleSpec("strain_proxy_target", "target", "none", target),
        SampleSpec("time_shuffle_control", "negative_control", "global_time_shuffle", time_shuffle(target, rng)),
        SampleSpec("phase_scramble_control", "negative_control", "fft_phase_scramble", phase_scramble(target, rng)),
        SampleSpec(
            "amplitude_preserving_surrogate_control",
            "negative_control",
            "iaaft_amplitude_spectrum_surrogate",
            amplitude_preserving_surrogate(target, rng, iterations=28),
        ),
        SampleSpec(
            "sliding_window_scramble_control",
            "negative_control",
            "event_window_sliding_scramble",
            sliding_window_scramble(record, config, target, rng),
        ),
        SampleSpec(
            "partial_overlap_scramble_control",
            "negative_control",
            "partial_event_overlap_scramble",
            partial_overlap_scramble(record, config, target, rng),
        ),
        SampleSpec(
            "chirp_reversal_control",
            "negative_control",
            "event_window_chirp_reversal",
            chirp_reversal(record, config, target),
        ),
        SampleSpec(
            "envelope_locked_phase_control",
            "negative_control",
            "event_envelope_preserved_phase_scramble",
            envelope_locked_phase_control(record, config, target, rng),
        ),
        SampleSpec(
            "multi_event_timing_randomization_control",
            "negative_control",
            "multi_event_timing_randomization",
            multi_event_timing_randomization(record, config, target, rng),
        ),
        SampleSpec(
            "micro_chunk_scramble_control",
            "negative_control",
            "event_window_micro_chunk_scramble",
            micro_chunk_scramble(record, config, target, rng),
        ),
    ]


def event_indices(record: StrainRecord, config: ProxyConfig) -> np.ndarray:
    return np.flatnonzero(event_masks(record.time, record.event_time, config)["event"])


def sliding_window_scramble(
    record: StrainRecord,
    config: ProxyConfig,
    strain: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.array(strain, copy=True)
    indices = event_indices(record, config)
    if indices.size < 32:
        return time_shuffle(out, rng)
    window = max(8, int(round(0.055 * record.sample_rate)))
    hop = max(4, window // 3)
    starts = list(range(0, max(1, indices.size - window), hop))
    if len(starts) < 2:
        out[indices] = out[indices[::-1]]
        return normalize_signal(out)
    pieces = [out[indices[start : start + window]].copy() for start in starts]
    order = rng.permutation(len(pieces))
    rebuilt = np.array(out[indices], copy=True)
    weights = np.zeros_like(rebuilt)
    acc = np.zeros_like(rebuilt)
    cursor = 0
    for piece_idx in order:
        start = starts[cursor % len(starts)]
        piece = pieces[int(piece_idx)]
        if rng.random() < 0.5:
            piece = piece[::-1]
        stop = min(start + piece.size, rebuilt.size)
        acc[start:stop] += piece[: stop - start]
        weights[start:stop] += 1.0
        cursor += 1
    keep = weights > 0
    rebuilt[keep] = acc[keep] / weights[keep]
    out[indices] = rebuilt
    return normalize_signal(out)


def partial_overlap_scramble(
    record: StrainRecord,
    config: ProxyConfig,
    strain: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.array(strain, copy=True)
    indices = event_indices(record, config)
    if indices.size < 24:
        return time_shuffle(out, rng)
    thirds = np.array_split(indices, 3)
    early, middle, late = thirds
    # Keep the middle energetic overlap but swap/reverse the flanks.
    flank_size = min(early.size, late.size)
    out[early[:flank_size]] = strain[late[:flank_size]][::-1]
    out[late[:flank_size]] = strain[early[:flank_size]][::-1]
    middle_values = np.array(strain[middle], copy=True)
    jitter = rng.normal(scale=0.08 * max(float(np.std(middle_values)), 1.0e-12), size=middle_values.size)
    out[middle] = middle_values + jitter
    return normalize_signal(out)


def chirp_reversal(record: StrainRecord, config: ProxyConfig, strain: np.ndarray) -> np.ndarray:
    out = np.array(strain, copy=True)
    indices = event_indices(record, config)
    if indices.size < 16:
        return normalize_signal(out[::-1])
    out[indices] = out[indices][::-1]
    return normalize_signal(out)


def envelope_locked_phase_control(
    record: StrainRecord,
    config: ProxyConfig,
    strain: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    products = compute_products(record, config, strain)
    carrier = phase_scramble(products.filtered, rng)
    envelope = products.envelope / max(float(np.max(products.envelope)), 1.0e-12)
    baseline = np.array(strain, copy=True)
    indices = event_indices(record, config)
    if indices.size < 16:
        return phase_scramble(strain, rng)
    pre_mask = event_masks(record.time, record.event_time, config)["pre"]
    noise_scale = float(np.std(strain[pre_mask])) if np.any(pre_mask) else float(np.std(strain))
    baseline = rng.normal(scale=max(noise_scale, 1.0e-6), size=strain.size)
    baseline[indices] += 1.18 * envelope[indices] * carrier[indices]
    return normalize_signal(baseline)


def multi_event_timing_randomization(
    record: StrainRecord,
    config: ProxyConfig,
    strain: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.array(strain, copy=True)
    indices = event_indices(record, config)
    if indices.size < 16:
        return time_shuffle(out, rng)
    event_segment = np.array(strain[indices], copy=True)
    out[indices] = rng.normal(scale=max(float(np.std(strain)), 1.0e-12), size=indices.size)
    possible_centers = [
        record.event_time - 0.46,
        record.event_time - 0.30,
        record.event_time + 0.27,
        record.event_time + 0.43,
    ]
    rng.shuffle(possible_centers)
    for weight, center in zip((0.72, 0.52), possible_centers[:2], strict=True):
        target_indices = np.flatnonzero(np.abs(record.time - center) <= config.event_half_width * 0.58)
        if target_indices.size < 8:
            continue
        size = min(target_indices.size, event_segment.size)
        start = max(0, (event_segment.size - size) // 2)
        out[target_indices[:size]] += weight * event_segment[start : start + size]
    return normalize_signal(out)


def micro_chunk_scramble(
    record: StrainRecord,
    config: ProxyConfig,
    strain: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.array(strain, copy=True)
    indices = event_indices(record, config)
    if indices.size < 24:
        return time_shuffle(out, rng)
    chunk = max(3, int(round(0.010 * record.sample_rate)))
    chunks = [indices[start : start + chunk] for start in range(0, indices.size, chunk)]
    values = [np.array(out[item], copy=True) for item in chunks]
    rng.shuffle(values)
    rebuilt = np.concatenate([item[::-1] if rng.random() < 0.5 else item for item in values])
    out[indices] = rebuilt[: indices.size]
    return normalize_signal(out)


def summarize_event(
    event_id: str,
    record: StrainRecord,
    seed: int,
    snr: float,
    summaries: list[SampleSummary],
    config: ProxyConfig,
) -> dict[str, Any]:
    target = next(row for row in summaries if row.sample_role == "target")
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    control_scores = np.array([row.strain_proxy_score for row in controls], dtype=float)
    threshold = target.strain_proxy_score - config.pass_min_margin
    best = max(controls, key=lambda row: row.strain_proxy_score)
    event_window_controls = [
        row
        for row in controls
        if "event" in row.control_type or "chirp" in row.control_type or "timing" in row.control_type
    ]
    event_best = max(event_window_controls, key=lambda row: row.strain_proxy_score)
    return {
        "event_id": event_id,
        "source_label": record.source_label,
        "source_kind": record.source_kind,
        "seed": int(seed),
        "signal_to_noise": float(snr),
        "target_score": float(target.strain_proxy_score),
        "best_control_sample": best.sample,
        "best_control_type": best.control_type,
        "best_control_score": float(best.strain_proxy_score),
        "margin_over_best_control": float(target.strain_proxy_score - best.strain_proxy_score),
        "event_window_control_max": float(event_best.strain_proxy_score),
        "event_window_margin": float(target.strain_proxy_score - event_best.strain_proxy_score),
        "separation_z": float(
            (target.strain_proxy_score - float(np.mean(control_scores)))
            / max(float(np.std(control_scores)), 1.0e-9)
        ),
        "control_false_positive_rate": float(np.mean(control_scores >= threshold)),
    }


def run_hardening(config: HardeningConfig, args: argparse.Namespace) -> tuple[list[EventRun], dict[str, Any]]:
    runs: list[EventRun] = []
    for event_id, seed, snr, record in build_records(args, config):
        pconfig = proxy_config(config, seed, snr)
        samples = build_hardening_samples(record, pconfig, seed)
        summaries = [evaluate_sample(record, pconfig, sample) for sample in samples]
        summary = summarize_event(event_id, record, seed, snr, summaries, pconfig)
        runs.append(EventRun(event_id, seed, snr, record, summaries, summary))
    status = build_status(config, runs)
    return runs, status


def build_status(config: HardeningConfig, runs: list[EventRun]) -> dict[str, Any]:
    summaries = [run.summary for run in runs]
    margins = np.array([row["margin_over_best_control"] for row in summaries], dtype=float)
    event_window_margins = np.array([row["event_window_margin"] for row in summaries], dtype=float)
    target_scores = np.array([row["target_score"] for row in summaries], dtype=float)
    separation_z = np.array([row["separation_z"] for row in summaries], dtype=float)
    false_positive_rates = np.array([row["control_false_positive_rate"] for row in summaries], dtype=float)
    best_event = min(summaries, key=lambda row: row["margin_over_best_control"])
    best_control_event = max(summaries, key=lambda row: row["best_control_score"])

    aggregate_false_positive_rate = float(np.mean(false_positive_rates > 0.0))
    target_passes = bool(float(np.mean(target_scores)) >= config.pass_min_mean_target_score)
    controls_pass = bool(
        float(np.mean(margins)) >= config.pass_min_mean_margin
        and float(np.min(margins)) >= config.pass_min_min_margin
        and float(np.mean(event_window_margins)) >= config.pass_min_mean_event_window_margin
        and float(np.mean(separation_z)) >= config.pass_min_mean_separation_z
        and aggregate_false_positive_rate <= config.pass_max_control_false_positive_rate
    )
    if target_passes and controls_pass:
        bridge_status = "PASS"
        failure_mode = "NONE"
    elif target_passes and float(np.mean(separation_z)) >= 1.0:
        bridge_status = "MARGINAL"
        failure_mode = "EVENT_WINDOW_CONTROL_LEAK_REMAINS"
    else:
        bridge_status = "OPEN"
        failure_mode = "MULTI_EVENT_STRAIN_PROXY_NOT_STABLE"

    return {
        "phase": "65",
        "bridge_name": "gw_event_window_hardening",
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "success_definition": (
            "PASS means the Phase 64 strain-derived proxy survives a stricter event-window "
            "control ladder across multiple deterministic strain-like surrogate events or "
            "local strain files. It does not mean real gravitational memory, BMS charges, "
            "soft theorems, celestial holography, or S-matrix recovery."
        ),
        "config": asdict(config),
        "event_count": len(runs),
        "mean_target_score": float(np.mean(target_scores)),
        "std_target_score": float(np.std(target_scores)),
        "mean_margin_over_best_control": float(np.mean(margins)),
        "min_margin_over_best_control": float(np.min(margins)),
        "mean_event_window_margin": float(np.mean(event_window_margins)),
        "min_event_window_margin": float(np.min(event_window_margins)),
        "mean_separation_z": float(np.mean(separation_z)),
        "min_separation_z": float(np.min(separation_z)),
        "aggregate_control_false_positive_rate": aggregate_false_positive_rate,
        "weakest_event_id": str(best_event["event_id"]),
        "weakest_event_margin": float(best_event["margin_over_best_control"]),
        "highest_control_event_id": str(best_control_event["event_id"]),
        "highest_control_sample": str(best_control_event["best_control_sample"]),
        "highest_control_score": float(best_control_event["best_control_score"]),
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "66_real_gwosc_event_replicates_or_ward_identity_gate",
        "non_claims": [
            "no real gravitational-memory observable claim",
            "no BMS charge or soft-hair recovery claim",
            "no celestial holography recovery claim",
            "no real soft theorem claim",
            "no S-matrix reconstruction claim",
            "event-window hardened strain-derived proxy telemetry only",
        ],
    }


def write_diagnostics(path: Path, runs: list[EventRun]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DIAGNOSTIC_FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for run in runs:
            for row in run.summaries:
                payload = {
                    "event_id": run.event_id,
                    "source_label": run.record.source_label,
                    "source_kind": run.record.source_kind,
                    "seed": run.seed,
                    "signal_to_noise": run.signal_to_noise,
                    "sample": row.sample,
                    "sample_role": row.sample_role,
                    "control_type": row.control_type,
                    "strain_proxy_score": row.strain_proxy_score,
                    "memory_proxy_recoverability": row.memory_proxy_recoverability,
                    "event_localization_score": row.event_localization_score,
                    "chirp_order_score": row.chirp_order_score,
                    "composition_alignment_score": row.composition_alignment_score,
                    "permanence_score": row.permanence_score,
                    "event_energy_fraction": row.event_energy_fraction,
                    "peak_time_error": row.peak_time_error,
                    "centroid_slope_hz": row.centroid_slope_hz,
                }
                writer.writerow(payload)


def write_summary(path: Path, runs: list[EventRun]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for run in runs:
            writer.writerow({key: run.summary[key] for key in SUMMARY_FIELDNAMES})


def write_status(path: Path, status: dict[str, Any]) -> None:
    path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")


def write_report(path: Path, runs: list[EventRun], status: dict[str, Any]) -> None:
    summary_rows = "\n".join(
        "| {event_id} | {target_score:.6f} | {best_control_sample} | {best_control_score:.6f} | "
        "{margin_over_best_control:.6f} | {event_window_margin:.6f} | {separation_z:.6f} |".format(
            **run.summary
        )
        for run in runs
    )
    report = f"""# Phase 65 GW Event-Window Hardening

## Scope

This report is generated by `experiments/physics_bridge/gw_memory_entry_gate/run_gw_event_window_hardening.py`.

Phase 65 attacks the main Phase 64 leak: event-window controls that preserve timing-envelope structure while breaking event-internal chirp and composition structure. The default run uses deterministic GW150914-like surrogate replicates. Optional local strain files can be supplied with `--input-files`.

This is not a gravitational-memory detection. It is not a BMS charge, soft theorem, celestial holography, or S-matrix result.

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- event count: `{status['event_count']}`
- mean target score: `{status['mean_target_score']:.6f}`
- mean margin over best control: `{status['mean_margin_over_best_control']:.6f}`
- minimum margin over best control: `{status['min_margin_over_best_control']:.6f}`
- mean event-window margin: `{status['mean_event_window_margin']:.6f}`
- mean separation z: `{status['mean_separation_z']:.6f}`
- aggregate false-positive rate: `{status['aggregate_control_false_positive_rate']:.6f}`
- weakest event: `{status['weakest_event_id']}` at margin `{status['weakest_event_margin']:.6f}`
- highest control: `{status['highest_control_sample']}` in `{status['highest_control_event_id']}` at `{status['highest_control_score']:.6f}`

## Event Summary

| Event | Target | Best control | Control score | Margin | Event-window margin | z |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
{summary_rows}

## Control Ladder

- global time shuffle
- FFT phase scramble
- amplitude-preserving spectrum surrogate
- sliding event-window scramble
- partial event-overlap scramble
- chirp-order reversal inside the event window
- event-envelope-preserved phase scramble
- multi-event timing randomization
- micro-chunk event-window scramble

## Interpretation

The target rows are allowed to mean only this: the strain-derived memory/composition proxy remained more specific than the event-window hardening controls under the tested replicate set. This still does not provide an external gravitational-memory observable.

## Claim Gate

Allowed language:

- event-window hardened strain-derived proxy
- multi-event surrogate robustness
- structured strain-like event separated from stricter controls

Disallowed language:

- real gravitational memory detected
- BMS charge or soft hair recovered
- celestial holography recovered
- soft theorem recovered
- S-matrix reconstructed

## Generated Artifacts

- diagnostics CSV: `experiments/physics_bridge/gw_memory_entry_gate/event_window_hardening_diagnostics.csv`
- summary CSV: `experiments/physics_bridge/gw_memory_entry_gate/event_window_hardening_summary.csv`
- JSON: `experiments/physics_bridge/gw_memory_entry_gate/event_window_hardening_status.json`
- figures: `experiments/physics_bridge/gw_memory_entry_gate/figures/phase65_*.png`
"""
    path.write_text(report, encoding="utf-8")


def write_figures(output_dir: Path, runs: list[EventRun], status: dict[str, Any]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        (output_dir / "phase65_figures_not_generated.txt").write_text(
            "Matplotlib is not installed; rerun with `uv run --with numpy --with matplotlib ...`.\n",
            encoding="utf-8",
        )
        return

    figure_dir = output_dir / FIGURE_DIR.name
    figure_dir.mkdir(parents=True, exist_ok=True)
    event_ids = [run.event_id for run in runs]
    target_scores = np.array([run.summary["target_score"] for run in runs], dtype=float)
    control_scores = np.array([run.summary["best_control_score"] for run in runs], dtype=float)
    margins = np.array([run.summary["margin_over_best_control"] for run in runs], dtype=float)
    event_margins = np.array([run.summary["event_window_margin"] for run in runs], dtype=float)

    x = np.arange(len(runs))
    plt.figure(figsize=(10.5, 4.8))
    plt.plot(x, target_scores, marker="o", color="#1f7a5c", label="target")
    plt.plot(x, control_scores, marker="o", color="#a23e48", label="best control")
    plt.xticks(x, event_ids, rotation=35, ha="right", fontsize=8)
    plt.ylabel("strain proxy score")
    plt.title("Phase 65 target vs strongest control by event")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figure_dir / "phase65_event_scores.png", dpi=180)
    plt.close()

    plt.figure(figsize=(9.0, 4.8))
    plt.bar(x - 0.18, margins, width=0.36, color="#335c67", label="best-control margin")
    plt.bar(x + 0.18, event_margins, width=0.36, color="#e09f3e", label="event-window margin")
    plt.axhline(0.25, color="#9e2a2b", ls="--", lw=1.0, label="target mean bar")
    plt.xticks(x, event_ids, rotation=35, ha="right", fontsize=8)
    plt.ylabel("margin")
    plt.title("Phase 65 margin distribution")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figure_dir / "phase65_margin_distribution.png", dpi=180)
    plt.close()

    control_types: dict[str, list[float]] = {}
    for run in runs:
        target = next(row for row in run.summaries if row.sample_role == "target")
        for row in run.summaries:
            if row.sample_role != "negative_control":
                continue
            control_types.setdefault(row.sample, []).append(target.strain_proxy_score - row.strain_proxy_score)
    labels = list(control_types)
    means = [float(np.mean(control_types[label])) for label in labels]
    plt.figure(figsize=(11.0, 5.2))
    plt.bar(labels, means, color="#5b6770")
    plt.xticks(rotation=35, ha="right", fontsize=8)
    plt.ylabel("mean target-control margin")
    plt.title("Phase 65 control ladder erosion")
    plt.tight_layout()
    plt.savefig(figure_dir / "phase65_control_ladder.png", dpi=180)
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

    runs, status = run_hardening(config, args)
    write_diagnostics(output_dir / DIAGNOSTICS_PATH.name, runs)
    write_summary(output_dir / SUMMARY_PATH.name, runs)
    write_status(output_dir / STATUS_PATH.name, status)
    write_report(output_dir / REPORT_PATH.name, runs, status)
    write_figures(output_dir, runs, status)

    print("Phase 65 GW event-window hardening generated.")
    print(f"Diagnostics: {relative(output_dir / DIAGNOSTICS_PATH.name)}")
    print(f"Summary: {relative(output_dir / SUMMARY_PATH.name)}")
    print(f"JSON: {relative(output_dir / STATUS_PATH.name)}")
    print(f"Report: {relative(output_dir / REPORT_PATH.name)}")
    print(f"bridge_status: {status['bridge_status']} ({status['failure_mode']})")
    print(
        f"mean_target={status['mean_target_score']:.6f}; "
        f"mean_margin={status['mean_margin_over_best_control']:.6f}; "
        f"min_margin={status['min_margin_over_best_control']:.6f}; "
        f"mean_z={status['mean_separation_z']:.6f}"
    )


if __name__ == "__main__":
    main()
