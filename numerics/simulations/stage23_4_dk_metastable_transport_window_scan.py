#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import (
    CSV_FIELDS as BASE_FIELDS,
    WORK_PLOT_DIR,
    basin_id,
    field_grid_2d,
    load_runsheet,
    longest_constant_run,
    mean_pair_separation_series,
    ordered_peaks,
    plot_grades,
    plot_trace,
    plot_trajectories,
    selected_runs,
    weighted_center,
    crank_nicolson_states,
    packet_state,
)
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions, periodic_displacement

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_4_dk_metastable_transport_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_4_DK_Metastable_Transport_Window_Scan_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'triplet',
    'variant_label',
    'separation',
    'sigma_a',
    'sigma_b',
    'phase_offset',
    'grade_coupling_modulation',
    'kernel_stretch_x',
    'kernel_stretch_y',
    'persistence_gain',
    'grade_exchange_coherence',
    'transport_span',
    'post_collision_energy_asymmetry',
    'recurrence_indicator',
    'bounded_oscillatory_separation',
    'derived_label',
    'stable_transport_family',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.4 DK metastable transport window scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def pair_phase_difference(state: np.ndarray, packet_states: list[np.ndarray]) -> float:
    refs = []
    for packet in packet_states:
        overlap = np.vdot(packet, state)
        refs.append(float(np.angle(overlap)))
    if len(refs) < 2:
        return 0.0
    diff = refs[1] - refs[0]
    return float((diff + math.pi) % (2.0 * math.pi) - math.pi)


def grade_weights(vec: np.ndarray, block_sizes: tuple[int, ...]) -> list[float]:
    weights: list[float] = []
    start = 0
    total = float(np.sum(np.abs(vec) ** 2)) or 1.0
    for size in block_sizes:
        part = vec[start:start + size]
        weights.append(float(np.sum(np.abs(part) ** 2) / total))
        start += size
    return weights


def local_grade_mix(vec: np.ndarray, block_sizes: tuple[int, ...]) -> np.ndarray:
    n0, n1, _ = block_sizes
    q0 = vec[:n0]
    q1 = vec[n0:n0 + n1]
    out = np.zeros_like(vec)
    m = min(n0, n1)
    out[:m] += q1[:m]
    out[n0:n0 + m] += q0[:m]
    return out


def bounded_oscillation(series: np.ndarray, close_threshold: float) -> bool:
    if len(series) < 6:
        return False
    tail = series[len(series) // 2:]
    if tail.size == 0:
        return False
    below = tail <= close_threshold
    if np.mean(below) < 0.5:
        return False
    span = float(np.max(tail) - np.min(tail))
    return 0.02 <= span <= 0.18


def transport_span_series(center_series: np.ndarray) -> np.ndarray:
    if len(center_series) == 0:
        return np.zeros(0, dtype=float)
    cumulative = [0.0]
    for prev, nxt in zip(center_series[:-1], center_series[1:]):
        step = periodic_displacement(nxt[None, :], prev)[0]
        cumulative.append(cumulative[-1] + float(np.linalg.norm(step)))
    return np.asarray(cumulative, dtype=float)


def recurrence_indicator(series: np.ndarray, close_threshold: float) -> float:
    if len(series) < 6:
        return 0.0
    tail = series[len(series) // 3:]
    if tail.size < 4:
        return 0.0
    deriv = np.diff(tail)
    signs = np.sign(deriv)
    signs = signs[signs != 0.0]
    if signs.size < 2:
        return 0.0
    turns = np.sum(signs[1:] * signs[:-1] < 0.0)
    closeness = float(np.mean(tail <= close_threshold))
    return float((turns / max(1, signs.size - 1)) * closeness)


def grade_exchange_coherence(grade_hist: np.ndarray) -> tuple[float, np.ndarray]:
    if grade_hist.size == 0:
        return 0.0, np.zeros(0, dtype=float)
    signal = np.linalg.norm(grade_hist - grade_hist[0][None, :], axis=1)
    peak = float(np.max(signal))
    if peak <= 1.0e-12:
        return 0.0, signal
    coherence = float(np.mean(signal) / peak)
    return max(0.0, min(1.0, coherence)), signal


def final_peak_asymmetry(grid: np.ndarray) -> float:
    flat = np.sort(grid.ravel())[::-1]
    if flat.size < 2:
        return 0.0
    a = float(flat[0])
    b = float(flat[1])
    denom = a + b
    if denom <= 1.0e-12:
        return 0.0
    return float(abs(a - b) / denom)


def evolve_probe(
    base_operator,
    block_sizes: tuple[int, ...],
    packet_states: list[np.ndarray],
    grade_coupling_modulation: float,
    dt: float,
    steps: int,
) -> tuple[list[np.ndarray], list[float]]:
    state = np.sum(packet_states, axis=0)
    states = [state.copy()]
    phase_diffs = [pair_phase_difference(state, packet_states)]
    for _ in range(steps):
        acc = -(base_operator @ state)
        if grade_coupling_modulation > 0.0:
            acc += grade_coupling_modulation * local_grade_mix(state, block_sizes)
        state = state + dt * acc
        norm = np.linalg.norm(state)
        if norm > 1.0e-12:
            state = state / norm * np.linalg.norm(states[0])
        states.append(state.copy())
        phase_diffs.append(pair_phase_difference(state, packet_states))
    return states, phase_diffs


def classify_label(
    persistence_gain: float,
    transport_span: float,
    coherence: float,
    recurrence: float,
    bounded_sep: bool,
) -> str:
    if recurrence >= 0.18 and bounded_sep:
        return 'candidate proto-bound composite'
    if transport_span >= 0.18 and coherence >= 0.55 and persistence_gain <= 0.03:
        return 'metastable transport'
    if transport_span <= 0.08 and coherence < 0.35:
        return 'dispersive washout'
    return 'localized encounter'


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['resolution'])
    epsilon = 0.2
    amplitude = float(common['base_amplitude_scale']) * 0.5
    t_final = 0.9
    kick_cycles = 1.0

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes

    half_sep = 0.5 * float(run['separation'])
    centers = [
        np.array([0.50 - half_sep, 0.50], dtype=float),
        np.array([0.50 + half_sep, 0.50], dtype=float),
    ]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]
    sigmas = [float(run['sigma_a']), float(run['sigma_b'])]

    packet_states = [
        packet_state(
            positions=positions,
            block_sizes=block_sizes,
            center=centers[idx],
            sigma=sigmas[idx],
            amplitude=amplitude,
            phase_offset=0.0 if idx == 0 else float(run['phase_offset']),
            kick_vector=kicks[idx],
            kick_cycles=kick_cycles,
            active_grades=(0, 1),
        )
        for idx in range(2)
    ]

    dt = first_order_dt(complex_data.dirac_kahler, 0.35)
    steps = max(2, int(math.ceil(t_final / dt)))
    states, phase_diffs = evolve_probe(
        base_operator=complex_data.dirac_kahler,
        block_sizes=block_sizes,
        packet_states=packet_states,
        grade_coupling_modulation=float(run['grade_coupling_modulation']),
        dt=dt,
        steps=steps,
    )
    times = [idx * dt for idx in range(len(states))]

    grade_hist = np.asarray([grade_weights(state, block_sizes) for state in states], dtype=float)
    transfer_coherence, transfer_signal = grade_exchange_coherence(grade_hist)

    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(2)]
    raw_peak_counts: list[int] = []
    previous: list[np.ndarray] | None = None
    anchors = [center.copy() for center in centers]
    peak_value_series: list[tuple[float, float]] = []
    for idx, state in enumerate(states):
        grid = field_grid_2d(positions, state, resolution)
        current, raw_count = ordered_peaks(grid, 2, anchors if idx == 0 else None, previous)
        raw_peak_counts.append(raw_count)
        for track, point in zip(peak_tracks, current):
            track.append(point.copy())
        previous = [point.copy() for point in current]
        flat = np.sort(grid.ravel())[::-1]
        peak_value_series.append((float(flat[0]), float(flat[1]) if flat.size > 1 else 0.0))
    center_histories = [np.asarray(track, dtype=float) for track in peak_tracks]

    total_centers = np.asarray([weighted_center(positions, np.abs(state) ** 2) for state in states], dtype=float)
    basin_ids = [basin_id(center) for center in total_centers]
    coarse_basin_persistence = longest_constant_run(basin_ids)
    transport_series = transport_span_series(total_centers)
    transport_span = float(transport_series[-1]) if transport_series.size else 0.0

    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    mean_pair_distances, min_pair_distances = mean_pair_separation_series(center_histories)
    close_threshold = float(sigmas[0] + sigmas[1])
    composite_lifetime = float(np.sum(mean_pair_distances <= close_threshold) * sample_dt)
    binding_persistence = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    encounter_dwell_time = composite_lifetime
    min_separation = float(np.min(min_pair_distances)) if min_pair_distances.size else 0.0
    final_mean_separation = float(mean_pair_distances[-1]) if mean_pair_distances.size else 0.0
    post_collision_trend = final_mean_separation - min_separation
    tail = mean_pair_distances[max(0, int(0.75 * len(mean_pair_distances))):]
    bound_tail_fraction = float(np.mean(tail <= close_threshold)) if tail.size else 0.0
    initial_peak_count = int(raw_peak_counts[0]) if raw_peak_counts else 0

    final_grid = field_grid_2d(positions, states[-1], resolution)
    phase_array = np.asarray(phase_diffs, dtype=float)
    bounded_sep = int(bounded_oscillation(mean_pair_distances, close_threshold))
    recurrence = recurrence_indicator(mean_pair_distances, close_threshold)
    peak_asymmetry = final_peak_asymmetry(final_grid)
    grade_transfer_amplitude = float(np.max(transfer_signal)) if transfer_signal.size else 0.0

    row = {
        'run_id': run['run_id'],
        'geometry_id': 'tight_clustered_pair',
        'phase_id': run['variant_label'] if run['triplet'] == 'phase_transport_tuning' else 'clustered_family',
        'resolution': resolution,
        'graph_type': common['graph_type'],
        'packet_count': 2,
        'initial_peak_count': initial_peak_count,
        'collision_label': 'unresolved / mixed',
        'persistence_label': 'weak persistence gain' if encounter_dwell_time >= 0.15 * float(times[-1]) else 'no persistence gain',
        'composite_lifetime': composite_lifetime,
        'binding_persistence': binding_persistence,
        'corridor_dwell': 0.0,
        'coarse_basin_persistence': coarse_basin_persistence,
        'minimum_separation': min_separation,
        'final_mean_separation': final_mean_separation,
        'post_collision_separation_trend': post_collision_trend,
        'encounter_dwell_time': encounter_dwell_time,
        'deflection_angle_proxy': 0.0,
        'reflection_fraction': 0.0,
        'grade_transfer_amplitude': grade_transfer_amplitude,
        'omega0_weight_initial': float(grade_hist[0, 0]),
        'omega1_weight_initial': float(grade_hist[0, 1]),
        'omega2_weight_initial': float(grade_hist[0, 2]),
        'omega0_weight_final': float(grade_hist[-1, 0]),
        'omega1_weight_final': float(grade_hist[-1, 1]),
        'omega2_weight_final': float(grade_hist[-1, 2]),
        'localized_circulation_proxy': 0.0,
        'new_collision_class': 0,
        'gate_met': 0,
        'promoted_followup': 0,
        'notes': run['notes'],
        'triplet': run['triplet'],
        'variant_label': run['variant_label'],
        'separation': float(run['separation']),
        'sigma_a': float(run['sigma_a']),
        'sigma_b': float(run['sigma_b']),
        'phase_offset': float(run['phase_offset']),
        'grade_coupling_modulation': float(run['grade_coupling_modulation']),
        'kernel_stretch_x': float(run['kernel_stretch_x']),
        'kernel_stretch_y': float(run['kernel_stretch_y']),
        'persistence_gain': 0.0,
        'grade_exchange_coherence': transfer_coherence,
        'transport_span': transport_span,
        'post_collision_energy_asymmetry': peak_asymmetry,
        'recurrence_indicator': recurrence,
        'bounded_oscillatory_separation': bounded_sep,
        'derived_label': 'localized encounter',
        'stable_transport_family': 0,
    }

    detail = {
        'run': run,
        'metrics': row,
        'times': times,
        'mean_pair_distances': mean_pair_distances.tolist(),
        'minimum_pair_distances': min_pair_distances.tolist(),
        'phase_differences': phase_diffs,
        'grade_weights': grade_hist.tolist(),
        'peak_tracks': [hist.tolist() for hist in center_histories],
        'transport_span_series': transport_series.tolist(),
        'grade_exchange_signal': transfer_signal.tolist(),
        'peak_value_series': peak_value_series,
    }

    plot_paths: list[Path] = []
    traj_path = WORK_PLOT_DIR / f"stage23_4_{run['run_id']}_collision_trajectory.png"
    sep_path = WORK_PLOT_DIR / f"stage23_4_{run['run_id']}_separation_trace.png"
    grade_path = WORK_PLOT_DIR / f"stage23_4_{run['run_id']}_grade_weight_trace.png"
    span_path = WORK_PLOT_DIR / f"stage23_4_{run['run_id']}_transport_span_trace.png"
    plot_trajectories(traj_path, run['run_id'], center_histories, 'tight_clustered_pair')
    plot_trace(sep_path, run['run_id'], times, mean_pair_distances, 'mean peak separation', 'Stage 23.4 separation trace')
    plot_grades(grade_path, run['run_id'], times, grade_hist)
    plot_trace(span_path, run['run_id'], times, transport_series, 'transport span', 'Stage 23.4 transport-span trace')
    plot_paths.extend([traj_path, sep_path, grade_path, span_path])
    return row, detail, plot_paths


def plot_transport_panel(path: Path, details: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(12.0, 10.0), sharex=True, sharey=True)
    for ax, detail in zip(axes.flat, details):
        times = np.asarray(detail['times'], dtype=float)
        series = np.asarray(detail['transport_span_series'], dtype=float)
        ax.plot(times, series, color='tab:blue')
        ax.set_title(detail['run']['run_id'].replace('S23_4_', ''), fontsize=8)
        ax.grid(alpha=0.25)
    for ax in axes[-1, :]:
        ax.set_xlabel('time')
    for ax in axes[:, 0]:
        ax.set_ylabel('transport span')
    fig.suptitle('Stage 23.4 transport span vs time', fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_coherence_heatmap(path: Path, details: list[dict[str, Any]]) -> None:
    max_len = max(len(detail['grade_exchange_signal']) for detail in details)
    heat = np.zeros((len(details), max_len), dtype=float)
    for idx, detail in enumerate(details):
        signal = np.asarray(detail['grade_exchange_signal'], dtype=float)
        if signal.size == 0:
            continue
        peak = float(np.max(signal))
        normed = signal / peak if peak > 1.0e-12 else signal
        heat[idx, :normed.size] = normed
    labels = [detail['run']['run_id'].replace('S23_4_', '') for detail in details]
    fig, ax = plt.subplots(figsize=(11.0, 4.8))
    im = ax.imshow(heat, aspect='auto', cmap='magma', vmin=0.0, vmax=1.0)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel('time index')
    ax.set_title('Stage 23.4 grade coherence heatmap')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_classification_table(path: Path, rows: list[dict[str, Any]]) -> None:
    label_to_value = {
        'dispersive washout': 0,
        'localized encounter': 1,
        'metastable transport': 2,
        'candidate proto-bound composite': 3,
    }
    vals = np.asarray([[label_to_value[str(row['derived_label'])]] for row in rows], dtype=float)
    labels = [row['run_id'].replace('S23_4_', '') for row in rows]
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    im = ax.imshow(vals, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xticks([0], ['class'])
    ax.set_title('Stage 23.4 stability classification table')
    for idx, row in enumerate(rows):
        ax.text(0, idx, str(row['derived_label']).replace('candidate ', 'cand. '), ha='center', va='center', color='white', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], stable_triplets: list[str]) -> None:
    derived_counts = Counter(str(row['derived_label']) for row in rows)
    lines = [
        '# Stage 23.4 DK Metastable Transport Window Scan v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Total runs: {len(rows)}',
        '',
        '## Derived labels',
    ]
    for label, count in sorted(derived_counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend([
        '',
        '## Snapshot answers',
        f"- Does DK define a reproducible transport regime? {'Yes' if stable_triplets else 'No clear stable triplet yet.'}",
        f"- Is transport sensitive to phase geometry or width contrast? {'Yes' if any(row['triplet'] in {'width_contrast', 'phase_transport_tuning'} and row['derived_label'] == 'metastable transport' for row in rows) else 'Sensitivity stayed weak or mixed in this pass.'}",
        f"- Is there any signal of delayed persistence promotion? {'Only weak persistence-level effects.' if any(float(row['persistence_gain']) > 0.0 for row in rows) else 'No delayed persistence promotion signal appeared.'}",
        '',
        '## Stable metastable-transport triplets',
        f"- {stable_triplets if stable_triplets else 'none'}",
        '',
        '## Plots',
    ])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')
    NOTE_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    runs = selected_runs(runsheet['runs'], args.run_ids)
    common = runsheet['common_fields']

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        row, detail, run_plots = simulate_run(run, common)
        rows.append(row)
        details.append(detail)
        plot_paths.extend(run_plots)

    rows.sort(key=lambda item: item['run_id'])
    details.sort(key=lambda item: item['run']['run_id'])

    baseline = next((row for row in rows if row['run_id'] == 'S23_4_sep_tight_baseline'), rows[0])
    baseline_lifetime = float(baseline['composite_lifetime'])
    for row in rows:
        gain = float(row['composite_lifetime']) - baseline_lifetime
        row['persistence_gain'] = gain
        row['derived_label'] = classify_label(
            persistence_gain=gain,
            transport_span=float(row['transport_span']),
            coherence=float(row['grade_exchange_coherence']),
            recurrence=float(row['recurrence_indicator']),
            bounded_sep=bool(int(row['bounded_oscillatory_separation'])),
        )
        row['collision_label'] = row['derived_label']
        row['new_collision_class'] = int(row['derived_label'] in {'metastable transport', 'candidate proto-bound composite'})

    by_triplet: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_triplet[str(row['triplet'])].append(row)
    stable_triplets: list[str] = []
    for triplet, group in by_triplet.items():
        labels = Counter(str(row['derived_label']) for row in group)
        if labels.get('metastable transport', 0) >= 2:
            stable_triplets.append(triplet)
            for row in group:
                if row['derived_label'] == 'metastable transport':
                    row['stable_transport_family'] = 1

    transport_panel = WORK_PLOT_DIR / 'stage23_4_transport_span_panel.png'
    coherence_map = WORK_PLOT_DIR / 'stage23_4_grade_coherence_heatmap.png'
    class_table = WORK_PLOT_DIR / 'stage23_4_stability_classification_table.png'
    plot_transport_panel(transport_panel, details)
    plot_coherence_heatmap(coherence_map, details)
    plot_classification_table(class_table, rows)
    plot_paths.extend([transport_panel, coherence_map, class_table])

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': {
            'derived_labels': dict(Counter(row['derived_label'] for row in rows)),
            'stable_triplets': stable_triplets,
        },
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage23_4_dk_metastable_transport',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, stable_triplets)


if __name__ == '__main__':
    main()
