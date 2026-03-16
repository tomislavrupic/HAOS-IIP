#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import (
    CSV_FIELDS as BASE_FIELDS,
    WORK_PLOT_DIR,
    field_grid_2d,
    mean_pair_separation_series,
    ordered_peaks,
    plot_trace,
)
from stage23_6_dk_phase_sensitivity_braid_exchange import (
    edge_grid,
    evolve,
    grade_exchange_signal,
    grade_weights,
    pair_phase_difference,
    packet_state_with_skew,
)
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_12_harmonic_detuning_continuum_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_12_Harmonic_Detuning_Continuum_Scan_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'run_order',
    'variant_label',
    'delta',
    'right_address',
    'topology_class',
    'flow_concentration_index',
    'address_selectivity_index',
    'address_sorting_score',
    'topology_survival_time',
    'transfer_asymmetry',
    'transfer_smear_index',
    'grade_exchange_coherence',
    'channel_count',
    'loop_count',
    'phase_alignment_metric',
    'continuity_band_label',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.12 harmonic detuning continuum scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_runs(runs: list[dict[str, Any]], run_ids: list[str]) -> list[dict[str, Any]]:
    if not run_ids:
        return runs
    wanted = set(run_ids)
    return [run for run in runs if str(run['run_id']) in wanted]


def grade_index_map(block_sizes: tuple[int, ...]) -> np.ndarray:
    labels = np.zeros(sum(block_sizes), dtype=int)
    start = 0
    for grade, size in enumerate(block_sizes):
        labels[start:start + size] = grade
        start += size
    return labels


def continuous_mod_distance(a: np.ndarray, b: np.ndarray, modulus: float) -> np.ndarray:
    diff = np.mod(np.abs(a - b), modulus)
    return np.minimum(diff, modulus - diff)


def compatibility_weight(distance: np.ndarray, eta_match: float, eta_mismatch: float) -> np.ndarray:
    clipped = np.clip(np.asarray(distance, dtype=float), 0.0, 2.0)
    weights = np.empty_like(clipped, dtype=float)
    near_mask = clipped <= 1.0
    weights[near_mask] = 1.0 + (eta_match - 1.0) * clipped[near_mask]
    weights[~near_mask] = eta_match + (eta_mismatch - eta_match) * (clipped[~near_mask] - 1.0)
    return weights


def build_address_labels(
    packet_states: list[np.ndarray],
    block_sizes: tuple[int, ...],
    left_address: float,
    right_address: float,
    support_floor_fraction: float,
    dominant_ratio_threshold: float,
) -> np.ndarray:
    packet_power = np.stack([np.abs(packet) ** 2 for packet in packet_states], axis=1)
    total_power = np.sum(packet_power, axis=1)
    max_power = np.max(total_power) if total_power.size else 0.0
    active = total_power >= float(support_floor_fraction) * max(max_power, 1.0e-12)
    dominant = np.argmax(packet_power, axis=1)
    dominant_power = packet_power[np.arange(packet_power.shape[0]), dominant]
    dominant_ratio = dominant_power / np.maximum(total_power, 1.0e-12)
    mixed = active & (dominant_ratio < float(dominant_ratio_threshold))
    grade_ids = grade_index_map(block_sizes)

    labels = np.full(total_power.shape, float(left_address), dtype=float)
    labels[active & (dominant == 1)] = float(right_address)
    labels[mixed] = np.where(dominant[mixed] == 0, float(left_address), float(right_address))
    labels[~active] = np.where((grade_ids[~active] % 2) == 0, float(left_address), float(right_address))
    return labels


def address_weighted_operator(
    base_operator: sp.csr_matrix,
    labels: np.ndarray,
    modulus: float,
    eta_match: float,
    eta_mismatch: float,
) -> sp.csr_matrix:
    coo = base_operator.tocoo()
    distances = continuous_mod_distance(labels[coo.row], labels[coo.col], modulus)
    weights = compatibility_weight(distances, eta_match, eta_mismatch)
    weighted = sp.csr_matrix((coo.data * weights, (coo.row, coo.col)), shape=base_operator.shape)
    return weighted.tocsr()


def address_activity_statistics(
    states: list[np.ndarray],
    operator: sp.csr_matrix,
    labels: np.ndarray,
    modulus: float,
    eta_match: float,
    eta_mismatch: float,
) -> tuple[float, float, float]:
    coo = operator.tocoo()
    mask = coo.row < coo.col
    if not np.any(mask):
        return 0.0, 0.0, 0.0
    row = coo.row[mask]
    col = coo.col[mask]
    distances = continuous_mod_distance(labels[row], labels[col], modulus)
    compat = compatibility_weight(distances, eta_match, eta_mismatch)
    activity = np.zeros(row.shape[0], dtype=float)
    for state in states:
        activity += np.abs(state[row]) * np.abs(state[col])
    total = float(np.sum(activity))
    if total <= 1.0e-12:
        return 0.0, 0.0, 0.0
    same_mask = distances <= 1.0e-9
    cross_mask = ~same_mask
    same_activity = float(np.sum(activity[same_mask]) / total)
    cross_activity = float(np.sum(activity[cross_mask]) / total)
    weighted_mean = float(np.sum(compat * activity) / total)
    min_weight = float(compatibility_weight(np.asarray([2.0]), eta_match, eta_mismatch)[0])
    selectivity = float(np.clip((weighted_mean - min_weight) / max(1.0 - min_weight, 1.0e-12), 0.0, 1.0))
    sorting = float(np.clip(same_activity - cross_activity, -1.0, 1.0))
    return same_activity, sorting, selectivity


def phase_alignment_metric(phase_diffs: list[float]) -> float:
    if not phase_diffs:
        return 0.0
    mean_abs = float(np.mean(np.abs(np.asarray(phase_diffs, dtype=float))))
    return float(np.clip(1.0 - mean_abs / math.pi, 0.0, 1.0))


def transfer_smear_index(concentration: float, coherence: float) -> float:
    return float(np.clip(1.0 - 0.5 * (float(concentration) + float(coherence)), 0.0, 1.0))


def classify_topology(concentration: float, coherence: float, braid_flag: bool, selectivity: float) -> str:
    if braid_flag and concentration >= 0.89 and coherence >= 0.44 and selectivity >= 0.88:
        return 'braid_like_exchange'
    if concentration >= 0.86 and coherence >= 0.38 and selectivity >= 0.72:
        return 'transfer_smeared'
    if concentration >= 0.82 and coherence >= 0.30:
        return 'mixed'
    return 'address_induced_smear'


def classify_trace_state(concentration: float, transfer_norm: float, braid_flag: bool, selectivity: float) -> str:
    if braid_flag and concentration >= 0.885 and transfer_norm >= 0.55 and selectivity >= 0.86:
        return 'braid_like_exchange'
    if concentration >= 0.855 and transfer_norm >= 0.36 and selectivity >= 0.70:
        return 'transfer_smeared'
    if concentration >= 0.80 and transfer_norm >= 0.20:
        return 'mixed'
    return 'address_induced_smear'


def topology_trace(
    edge_grids: list[np.ndarray],
    braid_votes: list[float],
    transfer_signal: np.ndarray,
    selectivity: float,
) -> tuple[list[str], list[float]]:
    traces: list[str] = []
    encoded: list[float] = []
    peak_transfer = float(np.max(transfer_signal)) if transfer_signal.size else 0.0
    mapping = {
        'braid_like_exchange': 0.0,
        'transfer_smeared': 1.0,
        'mixed': 2.0,
        'address_induced_smear': 3.0,
    }
    for idx, grid in enumerate(edge_grids):
        total = float(np.sum(grid))
        flat = np.sort(grid.ravel())[::-1]
        top_k = max(1, int(math.ceil(0.1 * flat.size)))
        concentration = float(np.sum(flat[:top_k]) / max(total, 1.0e-12))
        braid_flag = abs(float(np.mean(braid_votes[: idx + 1])) - 0.5) > 0.25 if braid_votes[: idx + 1] else False
        transfer_norm = float(transfer_signal[idx] / peak_transfer) if peak_transfer > 1.0e-12 and idx < transfer_signal.size else 0.0
        label = classify_trace_state(concentration, transfer_norm, braid_flag, selectivity)
        traces.append(label)
        encoded.append(mapping[label])
    return traces, encoded


def transfer_asymmetry(avg_edge_grid: np.ndarray) -> float:
    half = avg_edge_grid.shape[0] // 2
    left = float(np.sum(avg_edge_grid[:half, :]))
    right = float(np.sum(avg_edge_grid[half:, :]))
    total = left + right
    if total <= 1.0e-12:
        return 0.0
    return float(abs(left - right) / total)


def continuity_band_label(delta: float) -> str:
    if delta <= 0.125:
        return 'exact_band'
    if delta <= 0.375:
        return 'mild_band'
    if delta <= 0.625:
        return 'mid_band'
    if delta <= 0.875:
        return 'strong_band'
    return 'max_band'


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['resolution'])
    epsilon = float(common['epsilon'])
    sigma = float(common['mean_width'])
    amplitude = float(common['amplitude'])
    t_final = float(common['t_final'])
    dt_scale = float(common['dt_scale'])
    kick_cycles = float(common['kick_cycles'])
    phase_offset = float(common['phase_offset_fraction_of_pi']) * math.pi
    skew = float(common['local_phase_skew_fraction'])
    beta = float(common['beta'])
    separation = float(common['separation'])
    modulus = float(common['address_modulus'])
    eta_match = float(common['eta_match'])
    eta_mismatch = float(common['eta_mismatch'])
    support_floor_fraction = float(common['support_floor_fraction'])
    dominant_ratio_threshold = float(common['dominant_ratio_threshold'])
    delta = float(run['delta'])

    left_address = 0.0
    right_address = 2.0 * delta

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    n0, n1, _ = block_sizes

    centers = [
        np.array([0.50 - 0.5 * separation, 0.50], dtype=float),
        np.array([0.50 + 0.5 * separation, 0.50], dtype=float),
    ]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]
    packet_states = [
        packet_state_with_skew(
            positions=positions,
            block_sizes=block_sizes,
            center=centers[idx],
            sigma=sigma,
            amplitude=amplitude,
            phase_offset=0.0 if idx == 0 else phase_offset,
            kick_vector=kicks[idx],
            kick_cycles=kick_cycles,
            phase_skew_fraction=0.0 if idx == 0 else skew,
        )
        for idx in range(2)
    ]

    labels = build_address_labels(
        packet_states=packet_states,
        block_sizes=block_sizes,
        left_address=left_address,
        right_address=right_address,
        support_floor_fraction=support_floor_fraction,
        dominant_ratio_threshold=dominant_ratio_threshold,
    )

    base_operator = complex_data.dirac_kahler
    effective_operator = address_weighted_operator(base_operator, labels, modulus, eta_match, eta_mismatch)
    dt = first_order_dt(base_operator, dt_scale)

    psi0 = np.sum(packet_states, axis=0)
    steps = max(2, int(math.ceil(t_final / dt)))
    states = evolve(effective_operator, block_sizes, psi0, dt, steps, beta)
    times = [idx * dt for idx in range(len(states))]

    grade_hist = np.asarray([grade_weights(state, block_sizes) for state in states], dtype=float)
    transfer_signal, coherence = grade_exchange_signal(grade_hist)

    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(2)]
    previous: list[np.ndarray] | None = None
    anchors = [center.copy() for center in centers]
    braid_votes = []
    edge_grids = []
    phase_diffs = []
    raw_peak_counts: list[int] = []
    for idx, state in enumerate(states):
        grid = field_grid_2d(positions, state, resolution)
        current, raw_count = ordered_peaks(grid, 2, anchors if idx == 0 else None, previous)
        raw_peak_counts.append(raw_count)
        for track, point in zip(peak_tracks, current):
            track.append(point.copy())
        if len(current) == 2:
            braid_votes.append(float(current[0][0] > current[1][0]))
        previous = [point.copy() for point in current]
        edge_values = np.abs(state[n0:n0 + n1])
        edge_grids.append(edge_grid(complex_data.edge_midpoints, edge_values, resolution))
        phase_diffs.append(pair_phase_difference(state, packet_states))
    center_histories = [np.asarray(track, dtype=float) for track in peak_tracks]
    avg_edge_grid = np.mean(np.asarray(edge_grids, dtype=float), axis=0)

    mean_pair_distances, min_pair_distances = mean_pair_separation_series(center_histories)
    close_threshold = 2.0 * sigma
    composite_lifetime = float(np.sum(mean_pair_distances <= close_threshold) * dt)
    binding_persistence = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    final_mean_separation = float(mean_pair_distances[-1]) if mean_pair_distances.size else 0.0
    min_separation = float(np.min(min_pair_distances)) if min_pair_distances.size else 0.0
    post_collision_trend = final_mean_separation - min_separation

    total = float(np.sum(avg_edge_grid))
    flat = np.sort(avg_edge_grid.ravel())[::-1]
    top_k = max(1, int(math.ceil(0.1 * flat.size)))
    concentration = float(np.sum(flat[:top_k]) / max(total, 1.0e-12))
    braid_flag = abs(float(np.mean(braid_votes)) - 0.5) > 0.25 if braid_votes else False
    same_activity, sorting_score, selectivity = address_activity_statistics(
        states=states,
        operator=effective_operator,
        labels=labels,
        modulus=modulus,
        eta_match=eta_match,
        eta_mismatch=eta_mismatch,
    )
    topology = classify_topology(concentration, coherence, braid_flag, selectivity)
    smear_index = transfer_smear_index(concentration, coherence)
    topo_trace_labels, topo_trace_values = topology_trace(edge_grids, braid_votes, transfer_signal, selectivity)
    topology_survival_time = float(sum(1 for label in topo_trace_labels if label == 'braid_like_exchange') * dt)
    asymmetry = transfer_asymmetry(avg_edge_grid)
    phase_align = phase_alignment_metric(phase_diffs)

    row = {
        'run_id': str(run['run_id']),
        'geometry_id': 'tight_clustered_pair',
        'phase_id': str(run['variant_label']),
        'resolution': resolution,
        'graph_type': common['graph_type'],
        'packet_count': 2,
        'initial_peak_count': int(raw_peak_counts[0]) if raw_peak_counts else 0,
        'collision_label': topology,
        'persistence_label': 'weak persistence gain' if composite_lifetime >= 0.15 * float(times[-1]) else 'no persistence gain',
        'composite_lifetime': composite_lifetime,
        'binding_persistence': binding_persistence,
        'corridor_dwell': 0.0,
        'coarse_basin_persistence': 0.0,
        'minimum_separation': min_separation,
        'final_mean_separation': final_mean_separation,
        'post_collision_separation_trend': post_collision_trend,
        'encounter_dwell_time': composite_lifetime,
        'deflection_angle_proxy': 0.0,
        'reflection_fraction': 0.0,
        'grade_transfer_amplitude': float(np.max(transfer_signal)) if transfer_signal.size else 0.0,
        'omega0_weight_initial': float(grade_hist[0, 0]),
        'omega1_weight_initial': float(grade_hist[0, 1]),
        'omega2_weight_initial': float(grade_hist[0, 2]),
        'omega0_weight_final': float(grade_hist[-1, 0]),
        'omega1_weight_final': float(grade_hist[-1, 1]),
        'omega2_weight_final': float(grade_hist[-1, 2]),
        'localized_circulation_proxy': 0.0,
        'new_collision_class': int(topology == 'braid_like_exchange'),
        'gate_met': 0,
        'promoted_followup': 0,
        'run_order': int(run['run_order']),
        'variant_label': str(run['variant_label']),
        'delta': delta,
        'right_address': right_address,
        'topology_class': topology,
        'flow_concentration_index': concentration,
        'address_selectivity_index': selectivity,
        'address_sorting_score': sorting_score,
        'topology_survival_time': topology_survival_time,
        'transfer_asymmetry': asymmetry,
        'transfer_smear_index': smear_index,
        'grade_exchange_coherence': coherence,
        'channel_count': 0,
        'loop_count': 0,
        'phase_alignment_metric': phase_align,
        'continuity_band_label': continuity_band_label(delta),
        'notes': str(run['notes']),
    }

    detail = {
        'run': run,
        'metrics': row,
        'times': times,
        'topology_trace': topo_trace_labels,
        'topology_trace_values': topo_trace_values,
        'phase_differences': phase_diffs,
        'grade_weights': grade_hist.tolist(),
        'grade_exchange_trace': transfer_signal.tolist(),
        'flow_concentration_trace': [float(np.sum(np.sort(g.ravel())[::-1][:max(1, int(math.ceil(0.1 * g.size)))]) / max(np.sum(g), 1.0e-12)) for g in edge_grids],
    }

    plot_paths: list[Path] = []
    topology_trace_path = WORK_PLOT_DIR / f"stage_c0_12_{run['run_id']}_topology_vs_time_trace.png"
    plot_trace(
        topology_trace_path,
        str(run['run_id']),
        times[: len(topo_trace_values)],
        np.asarray(topo_trace_values, dtype=float),
        'topology code',
        'Stage C0.12 topology-vs-time trace',
    )
    plot_paths.append(topology_trace_path)
    return row, detail, plot_paths


def plot_topology_phase_diagram(path: Path, rows: list[dict[str, Any]]) -> None:
    mapping = {
        'braid_like_exchange': 0,
        'transfer_smeared': 1,
        'mixed': 2,
        'address_induced_smear': 3,
    }
    ordered = sorted(rows, key=lambda item: float(item['delta']))
    xs = [float(row['delta']) for row in ordered]
    ys = [mapping[str(row['topology_class'])] for row in ordered]
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.plot(xs, ys, marker='o', color='tab:blue')
    ax.set_yticks([0, 1, 2, 3], ['braid', 'transfer', 'mixed', 'smear'])
    ax.set_xlabel('detuning delta')
    ax.set_title('Stage C0.12 topology phase diagram')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_selectivity_decay(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=lambda item: float(item['delta']))
    xs = [float(row['delta']) for row in ordered]
    ys = [float(row['address_selectivity_index']) for row in ordered]
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.plot(xs, ys, marker='o', color='tab:purple')
    ax.set_xlabel('detuning delta')
    ax.set_ylabel('address selectivity index')
    ax.set_title('Stage C0.12 selectivity decay curve')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_sorting_degradation(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=lambda item: float(item['delta']))
    xs = [float(row['delta']) for row in ordered]
    ys = [float(row['address_sorting_score']) for row in ordered]
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.plot(xs, ys, marker='o', color='tab:orange')
    ax.set_xlabel('detuning delta')
    ax.set_ylabel('address sorting score')
    ax.set_title('Stage C0.12 sorting degradation curve')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_braid_survival(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=lambda item: float(item['delta']))
    xs = [float(row['delta']) for row in ordered]
    ys = [float(row['topology_survival_time']) for row in ordered]
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.bar(xs, ys, width=0.08, color='tab:green')
    ax.set_xlabel('detuning delta')
    ax.set_ylabel('braid survival time')
    ax.set_title('Stage C0.12 braid survival panel')
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def continuity_assessment(rows: list[dict[str, Any]]) -> tuple[str, str]:
    ordered = sorted(rows, key=lambda item: float(item['delta']))
    classes = [str(row['topology_class']) for row in ordered]
    selectivity = [float(row['address_selectivity_index']) for row in ordered]
    sorting = [float(row['address_sorting_score']) for row in ordered]
    unique_classes = list(dict.fromkeys(classes))
    transitions = sum(1 for left, right in zip(classes, classes[1:]) if left != right)
    selectivity_monotone = all(left >= right - 1.0e-6 for left, right in zip(selectivity, selectivity[1:]))
    sorting_monotone = all(left >= right - 1.0e-6 for left, right in zip(sorting, sorting[1:]))

    if transitions == 1 and selectivity_monotone:
        threshold_idx = next(idx for idx, (left, right) in enumerate(zip(classes, classes[1:])) if left != right)
        threshold_delta = float(ordered[threshold_idx + 1]['delta'])
        if sorting_monotone:
            return 'sharp_transition_threshold', f'Topology stays in one protected class and then switches once near delta={threshold_delta:.3f}, while selectivity and sorting remain ordered across the ladder.'
        return 'sharp_transition_threshold', f'Topology stays in one protected class and then switches once near delta={threshold_delta:.3f}. Selectivity decays monotonically, while the sorting curve shows secondary texture inside the protected branch rather than a second topology transition.'
    if transitions >= 2 and selectivity_monotone and sorting_monotone:
        return 'multi_plateau_resonance_ladder', 'Topology changes in more than one stable band while selectivity and sorting decay remain ordered, indicating a multi-plateau resonance ladder.'
    if transitions <= 1 and selectivity_monotone and len(unique_classes) == 1:
        return 'smooth_continuum_decay', 'The topology label stays flat while selectivity and sorting decay smoothly, so the visible response is continuous degradation rather than a discrete switch.'
    return 'incoherent_response', 'Topology and detuning do not organize into a clean threshold or plateau pattern, so the response is best described as incoherent.'


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], assessment: str, note: str) -> None:
    lines = [
        '# Stage C0.12 Harmonic Detuning Continuum Scan v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Continuity assessment: `{assessment}`',
        '',
        'Per-run summary:',
        '',
    ]
    ordered = sorted(rows, key=lambda item: float(item['delta']))
    for row in ordered:
        lines.extend([
            f"- `{row['run_id']}`",
            f"  - delta: `{row['delta']:.3f}`",
            f"  - topology class: `{row['topology_class']}`",
            f"  - flow concentration index: `{row['flow_concentration_index']:.4f}`",
            f"  - address selectivity index: `{row['address_selectivity_index']:.4f}`",
            f"  - address sorting score: `{row['address_sorting_score']:.4f}`",
            f"  - topology survival time: `{row['topology_survival_time']:.4f}`",
            f"  - transfer asymmetry: `{row['transfer_asymmetry']:.4f}`",
        ])
    lines.extend([
        '',
        'Continuity note:',
        f'- {note}',
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

    rows.sort(key=lambda item: float(item['delta']))
    details.sort(key=lambda item: float(item['metrics']['delta']))

    phase_path = WORK_PLOT_DIR / 'stage_c0_12_topology_phase_diagram.png'
    selectivity_path = WORK_PLOT_DIR / 'stage_c0_12_selectivity_decay_curve.png'
    sorting_path = WORK_PLOT_DIR / 'stage_c0_12_sorting_degradation_curve.png'
    survival_path = WORK_PLOT_DIR / 'stage_c0_12_braid_survival_panel.png'
    plot_topology_phase_diagram(phase_path, rows)
    plot_selectivity_decay(selectivity_path, rows)
    plot_sorting_degradation(sorting_path, rows)
    plot_braid_survival(survival_path, rows)
    plot_paths.extend([phase_path, selectivity_path, sorting_path, survival_path])

    assessment, note = continuity_assessment(rows)
    summary = {
        'topology_counts': dict(Counter(str(row['topology_class']) for row in rows)),
        'assessment': assessment,
        'continuity_note': note,
    }
    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': summary,
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage_c0_12_harmonic_detuning_continuum_scan',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, assessment, note)


if __name__ == '__main__':
    main()
