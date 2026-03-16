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
    field_grid_2d,
    mean_pair_separation_series,
    ordered_peaks,
    weighted_center,
)
from stage23_6_dk_phase_sensitivity_braid_exchange import connected_components, edge_grid, grade_exchange_signal
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions, periodic_displacement

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_9_dk_collision_phenomenology_consolidation_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_9_DK_Collision_Phenomenology_Consolidation_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage23_9_consolidation'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'family_id',
    'family_label',
    'phase_detuning_fraction_of_pi',
    'kernel_radius_label',
    'kernel_radius_offset',
    'resolution',
    'refinement_label',
    'topology_class',
    'braid_like_exchange',
    'transfer_smeared',
    'localized_capture',
    'unresolved_mixed',
    'braid_survival_time',
    'flow_concentration_index',
    'grade_exchange_coherence',
    'recurrence_indicator',
    'midpoint_band_membership_flag',
    'metastable_transport_flag',
    'structural_irreversibility_score',
    'topology_return_error',
    'refinement_survival_ratio',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.9 DK collision phenomenology consolidation lattice.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def packet_state_weighted(
    positions: np.ndarray,
    block_sizes: tuple[int, ...],
    center: np.ndarray,
    sigma: float,
    amplitude: float,
    phase_offset: float,
    kick_vector: np.ndarray,
    kick_cycles: float,
    grade0_scale: float,
    grade1_scale: float,
) -> np.ndarray:
    parts: list[np.ndarray] = []
    start = 0
    for grade_idx, size in enumerate(block_sizes):
        coords = positions[start:start + size]
        delta = periodic_displacement(coords, center)
        r2 = np.sum(delta * delta, axis=1)
        profile = np.exp(-0.5 * r2 / max(sigma * sigma, 1.0e-12))
        phase = np.exp(1j * (2.0 * math.pi * kick_cycles * (delta @ kick_vector) + phase_offset))
        if grade_idx == 0:
            part = float(grade0_scale) * profile * phase
        elif grade_idx == 1:
            part = float(grade1_scale) * profile * phase
        else:
            part = np.zeros(size, dtype=complex)
        parts.append(np.asarray(part, dtype=complex))
        start += size
    psi = np.concatenate(parts)
    psi /= max(np.linalg.norm(psi), 1.0e-12)
    return float(amplitude) * psi


def local_grade_mix(vec: np.ndarray, block_sizes: tuple[int, ...]) -> np.ndarray:
    n0, n1, _ = block_sizes
    q0 = vec[:n0]
    q1 = vec[n0:n0 + n1]
    out = np.zeros_like(vec)
    m = min(n0, n1)
    out[:m] += q1[:m]
    out[n0:n0 + m] += q0[:m]
    return out


def evolve_signed(
    base_operator,
    block_sizes: tuple[int, ...],
    state0: np.ndarray,
    dt: float,
    steps: int,
    beta: float,
    direction: float,
) -> list[np.ndarray]:
    state = state0.copy()
    states = [state.copy()]
    base_norm = float(np.linalg.norm(state0))
    for _ in range(steps):
        acc = -direction * (base_operator @ state)
        if beta > 0.0:
            acc += direction * beta * local_grade_mix(state, block_sizes)
        state = state + dt * acc
        norm = float(np.linalg.norm(state))
        if not np.isfinite(norm) or norm <= 1.0e-12:
            raise RuntimeError('numerical instability indicator triggered')
        state = state / norm * base_norm
        states.append(state.copy())
    return states


def overlap_phase_alignment(state: np.ndarray, packet_states: list[np.ndarray]) -> float:
    phases: list[float] = []
    for packet in packet_states:
        phases.append(float(np.angle(np.vdot(packet, state))))
    if len(phases) < 2:
        return 0.0
    diffs = []
    anchor = phases[0]
    for phase in phases[1:]:
        diff = (phase - anchor + math.pi) % (2.0 * math.pi) - math.pi
        diffs.append(abs(diff))
    mean_abs = float(np.mean(diffs)) if diffs else 0.0
    return float(np.clip(1.0 - mean_abs / math.pi, 0.0, 1.0))


def flow_concentration(grid: np.ndarray) -> float:
    total = float(np.sum(grid))
    if total <= 1.0e-12:
        return 0.0
    flat = np.sort(grid.ravel())[::-1]
    top_k = max(1, int(math.ceil(0.1 * flat.size)))
    return float(np.sum(flat[:top_k]) / total)


def count_channels(grid: np.ndarray) -> int:
    threshold = max(0.15 * float(np.max(grid)), float(np.mean(grid) + np.std(grid)))
    mask = grid >= threshold
    return connected_components(mask)


def classify_topology(
    packet_count: int,
    braid_flag: bool,
    concentration: float,
    coherence: float,
    grade_amplitude: float,
    recurrence: float,
    close_fraction: float,
) -> str:
    if packet_count == 2 and braid_flag and concentration >= 0.84 and coherence >= 0.40:
        return 'braid_like_exchange'
    if packet_count == 2 and not braid_flag and 0.66 <= concentration <= 0.74 and coherence >= 0.40:
        return 'localized_capture'
    if close_fraction >= 0.58 and recurrence >= 0.08 and concentration >= 0.73:
        return 'localized_capture'
    if coherence >= 0.33 and grade_amplitude >= 0.09:
        return 'transfer_smeared'
    return 'unresolved_mixed'


def phase_offset_for_packet(packet_idx: int, packet_count: int, phase_offset: float) -> float:
    if packet_idx == 0:
        return 0.0
    if packet_count == 2:
        return phase_offset
    if packet_idx == 1:
        return phase_offset
    return -phase_offset


def analyze_states(
    states: list[np.ndarray],
    packet_states: list[np.ndarray],
    positions: np.ndarray,
    edge_midpoints: np.ndarray,
    block_sizes: tuple[int, ...],
    resolution: int,
    anchors: list[np.ndarray],
    sigma: float,
) -> dict[str, Any]:
    n0, n1, _ = block_sizes
    packet_count = len(packet_states)
    dt = 1.0
    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(packet_count)]
    previous: list[np.ndarray] | None = None
    edge_grids: list[np.ndarray] = []
    flow_trace: list[float] = []
    align_trace: list[float] = []
    raw_peak_counts: list[int] = []
    braid_votes: list[float] = []

    grade_hist = []
    for idx, state in enumerate(states):
        grid = field_grid_2d(positions, state, resolution)
        current, raw_count = ordered_peaks(grid, packet_count, anchors if idx == 0 else None, previous)
        raw_peak_counts.append(raw_count)
        for track, point in zip(peak_tracks, current):
            track.append(point.copy())
        previous = [point.copy() for point in current]
        if packet_count == 2 and len(current) == 2:
            braid_votes.append(float(current[0][0] > current[1][0]))
        edge_values = np.abs(state[n0:n0 + n1])
        egrid = edge_grid(edge_midpoints, edge_values, resolution)
        edge_grids.append(egrid)
        flow_trace.append(flow_concentration(egrid))
        align_trace.append(overlap_phase_alignment(state, packet_states))
        total = float(np.sum(np.abs(state) ** 2)) or 1.0
        grade_hist.append([
            float(np.sum(np.abs(state[:n0]) ** 2) / total),
            float(np.sum(np.abs(state[n0:n0 + n1]) ** 2) / total),
            float(np.sum(np.abs(state[n0 + n1:]) ** 2) / total),
        ])

    center_histories = [np.asarray(track, dtype=float) for track in peak_tracks]
    mean_pair_distances, min_pair_distances = mean_pair_separation_series(center_histories)
    close_threshold = 2.0 * sigma
    close_fraction = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0

    grade_hist_arr = np.asarray(grade_hist, dtype=float)
    transfer_signal, coherence = grade_exchange_signal(grade_hist_arr)
    grade_amplitude = float(np.max(transfer_signal)) if transfer_signal.size else 0.0
    recurrence = 0.0
    if mean_pair_distances.size >= 6:
        tail = mean_pair_distances[len(mean_pair_distances) // 3:]
        deriv = np.diff(tail)
        signs = np.sign(deriv)
        signs = signs[signs != 0.0]
        if signs.size >= 2:
            turns = np.sum(signs[1:] * signs[:-1] < 0.0)
            recurrence = float((turns / max(1, signs.size - 1)) * close_fraction)

    braid_flag = abs(float(np.mean(braid_votes)) - 0.5) > 0.25 if braid_votes else False
    avg_edge_grid = np.mean(np.asarray(edge_grids, dtype=float), axis=0)
    concentration = flow_concentration(avg_edge_grid)
    channels = count_channels(avg_edge_grid)

    topology_class = classify_topology(
        packet_count=packet_count,
        braid_flag=braid_flag,
        concentration=concentration,
        coherence=coherence,
        grade_amplitude=grade_amplitude,
        recurrence=recurrence,
        close_fraction=close_fraction,
    )

    instant_labels: list[str] = []
    max_grade = max(grade_amplitude, 1.0e-12)
    for idx, egrid in enumerate(edge_grids):
        local_conc = flow_concentration(egrid)
        local_sep = float(mean_pair_distances[idx]) if idx < len(mean_pair_distances) else float(mean_pair_distances[-1])
        local_close = 1.0 if local_sep <= close_threshold else 0.0
        local_coh = float(np.clip(float(transfer_signal[idx]) / max_grade, 0.0, 1.0)) if idx < len(transfer_signal) else coherence
        local_braid = braid_flag and packet_count == 2 and idx < len(braid_votes) and abs(braid_votes[idx] - 0.5) > 0.25
        instant_labels.append(
            classify_topology(
                packet_count=packet_count,
                braid_flag=local_braid,
                concentration=local_conc,
                coherence=local_coh,
                grade_amplitude=grade_amplitude,
                recurrence=recurrence,
                close_fraction=local_close,
            )
        )

    braid_survival_time = float(sum(label == 'braid_like_exchange' for label in instant_labels))
    return {
        'topology_class': topology_class,
        'flow_concentration_index': concentration,
        'grade_exchange_coherence': coherence,
        'recurrence_indicator': recurrence,
        'channel_count': channels,
        'braid_flag': braid_flag,
        'close_fraction': close_fraction,
        'braid_survival_time': braid_survival_time,
        'mean_pair_distances': mean_pair_distances.tolist(),
        'minimum_pair_distances': min_pair_distances.tolist(),
        'flow_trace': flow_trace,
        'alignment_trace': align_trace,
        'grade_exchange_trace': transfer_signal.tolist(),
        'topology_trace': instant_labels,
    }


def simulate_run(
    family: dict[str, Any],
    phase_fraction: float,
    radius_offset: int,
    resolution: int,
    common: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    epsilon = float(common['base_epsilon']) + float(radius_offset) / float(resolution)
    beta = float(common['beta'])
    phase_offset = float(phase_fraction) * math.pi
    kick_cycles = 1.0
    amplitude = 0.5

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    anchors = [np.asarray(center, dtype=float) for center in family['centers']]

    packet_states = []
    for idx, center in enumerate(family['centers']):
        packet_states.append(
            packet_state_weighted(
                positions=positions,
                block_sizes=block_sizes,
                center=np.asarray(center, dtype=float),
                sigma=float(family['sigma']),
                amplitude=amplitude,
                phase_offset=phase_offset_for_packet(idx, int(family['packet_count']), phase_offset),
                kick_vector=np.asarray(family['kicks'][idx], dtype=float),
                kick_cycles=kick_cycles,
                grade0_scale=float(family['grade0_scale']),
                grade1_scale=float(family['grade1_scale']),
            )
        )

    psi0 = np.sum(packet_states, axis=0)
    dt = first_order_dt(complex_data.dirac_kahler, float(common['dt_scale']))
    steps = max(2, int(math.ceil(float(common['t_final']) / dt)))

    forward_states = evolve_signed(complex_data.dirac_kahler, block_sizes, psi0, dt, steps, beta, direction=1.0)
    reverse_states = evolve_signed(complex_data.dirac_kahler, block_sizes, forward_states[-1], dt, steps, beta, direction=-1.0)

    forward_metrics = analyze_states(
        forward_states,
        packet_states,
        positions,
        complex_data.edge_midpoints,
        block_sizes,
        resolution,
        anchors,
        float(family['sigma']),
    )
    reverse_metrics = analyze_states(
        reverse_states,
        packet_states,
        positions,
        complex_data.edge_midpoints,
        block_sizes,
        resolution,
        anchors,
        float(family['sigma']),
    )

    ftrace = np.asarray(forward_metrics['flow_trace'], dtype=float)
    rtrace = np.asarray(reverse_metrics['flow_trace'], dtype=float)
    structural_irreversibility_score = float(np.mean(np.abs(ftrace - rtrace[::-1]))) if ftrace.size and rtrace.size else 0.0

    returned_state = reverse_states[-1]
    return_alignment = float(np.linalg.norm(returned_state - psi0) / max(np.linalg.norm(psi0), 1.0e-12))
    topology_return_error = int(reverse_metrics['topology_class'] != forward_metrics['topology_class'])
    midpoint_flag = int(
        family['family_id'] == 'clustered_dk_baseline'
        and 0.45 <= float(phase_fraction) <= 0.575
        and forward_metrics['topology_class'] == 'transfer_smeared'
    )
    metastable_transport_flag = int(forward_metrics['braid_survival_time'] * dt > 0.2 and forward_metrics['recurrence_indicator'] > 0.0)

    run_id = (
        f"S23_9_{family['family_id']}"
        f"_p{str(phase_fraction).replace('.', 'p')}"
        f"_k{radius_offset:+d}"
        f"_n{resolution}"
    )
    row = {
        'run_id': run_id,
        'family_id': family['family_id'],
        'family_label': family['family_label'],
        'phase_detuning_fraction_of_pi': float(phase_fraction),
        'kernel_radius_label': {-1: 'r_minus_1', 0: 'r', 1: 'r_plus_1'}[int(radius_offset)],
        'kernel_radius_offset': int(radius_offset),
        'resolution': int(resolution),
        'refinement_label': 'base' if int(resolution) == int(common['base_resolution']) else 'refined',
        'topology_class': forward_metrics['topology_class'],
        'braid_like_exchange': int(forward_metrics['topology_class'] == 'braid_like_exchange'),
        'transfer_smeared': int(forward_metrics['topology_class'] == 'transfer_smeared'),
        'localized_capture': int(forward_metrics['topology_class'] == 'localized_capture'),
        'unresolved_mixed': int(forward_metrics['topology_class'] == 'unresolved_mixed'),
        'braid_survival_time': float(forward_metrics['braid_survival_time'] * dt),
        'flow_concentration_index': float(forward_metrics['flow_concentration_index']),
        'grade_exchange_coherence': float(forward_metrics['grade_exchange_coherence']),
        'recurrence_indicator': float(forward_metrics['recurrence_indicator']),
        'midpoint_band_membership_flag': midpoint_flag,
        'metastable_transport_flag': metastable_transport_flag,
        'structural_irreversibility_score': structural_irreversibility_score,
        'topology_return_error': topology_return_error,
        'refinement_survival_ratio': 0.0,
    }
    detail = {
        'family': family,
        'phase_detuning_fraction_of_pi': float(phase_fraction),
        'kernel_radius_offset': int(radius_offset),
        'epsilon': epsilon,
        'dt': dt,
        'steps': steps,
        'return_alignment_error': return_alignment,
        'metrics': row,
        'forward': forward_metrics,
        'reverse': reverse_metrics,
    }

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f'{run_id}_trace.png'
    times = np.arange(len(forward_metrics['flow_trace']), dtype=float) * dt
    fig, axes = plt.subplots(3, 1, figsize=(6.6, 7.6), sharex=True)
    axes[0].plot(times, np.asarray(forward_metrics['flow_trace'], dtype=float), color='tab:blue')
    axes[0].set_ylabel('flow')
    axes[0].set_title(run_id)
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, np.asarray(forward_metrics['grade_exchange_trace'], dtype=float), color='tab:orange')
    axes[1].set_ylabel('grade exch.')
    axes[1].grid(alpha=0.25)
    topo_map = {'braid_like_exchange': 0, 'transfer_smeared': 1, 'localized_capture': 2, 'unresolved_mixed': 3}
    axes[2].plot(times, [topo_map[label] for label in forward_metrics['topology_trace']], color='tab:green')
    axes[2].set_ylabel('topology')
    axes[2].set_yticks([0, 1, 2, 3], ['braid', 'smeared', 'capture', 'mixed'])
    axes[2].set_xlabel('time')
    axes[2].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)
    return row, detail, plot_paths


def apply_refinement_survival(rows: list[dict[str, Any]]) -> None:
    grouped: dict[tuple[str, float, int], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = (str(row['family_id']), float(row['phase_detuning_fraction_of_pi']), int(row['kernel_radius_offset']))
        grouped[key].append(row)
    for group in grouped.values():
        by_res = {int(row['resolution']): row for row in group}
        if len(by_res) == 2:
            labels = {str(row['topology_class']) for row in by_res.values()}
            ratio = 1.0 if len(labels) == 1 else 0.0
            for row in group:
                row['refinement_survival_ratio'] = ratio


def classify_outcome(rows: list[dict[str, Any]]) -> str:
    clustered = [row for row in rows if row['family_id'] == 'clustered_dk_baseline']
    corridor = [row for row in rows if row['family_id'] == 'corridor_transport_family']
    triad = [row for row in rows if row['family_id'] == 'triad_competition_family']

    clustered_braid = sum(int(row['braid_like_exchange']) for row in clustered)
    clustered_midpoint_smear = sum(int(row['midpoint_band_membership_flag']) for row in clustered)
    transport_hits = sum(int(row['metastable_transport_flag']) for row in rows)
    corridor_capture = sum(int(row['localized_capture']) for row in corridor)
    triad_variety = len({str(row['topology_class']) for row in triad})
    irreversibility_pockets = sum(1 for row in rows if float(row['structural_irreversibility_score']) >= 0.01)

    if clustered_braid >= 12 and clustered_midpoint_smear >= 3 and transport_hits == 0 and irreversibility_pockets >= 1 and (corridor_capture >= 6 or triad_variety >= 2):
        return 'phenomenology_consolidated'
    if clustered_braid >= 4 and transport_hits == 0:
        return 'partial_structural_instability'
    return 'topology_family_inconsistent'


def plot_topology_stability_map(path: Path, rows: list[dict[str, Any]]) -> None:
    families = ['clustered_dk_baseline', 'corridor_transport_family', 'triad_competition_family']
    phases = sorted({float(row['phase_detuning_fraction_of_pi']) for row in rows})
    radii = [-1, 0, 1]
    value_map = {'absent': 0, 'fragile': 1, 'robust': 2}
    grid = np.zeros((len(families), len(phases) * len(radii)), dtype=float)
    for fi, family in enumerate(families):
        for pi, phase in enumerate(phases):
            for ri, radius in enumerate(radii):
                subset = [
                    row for row in rows
                    if row['family_id'] == family
                    and float(row['phase_detuning_fraction_of_pi']) == phase
                    and int(row['kernel_radius_offset']) == radius
                ]
                braid_count = sum(int(row['braid_like_exchange']) for row in subset)
                if braid_count == 2:
                    state = 'robust'
                elif braid_count == 1:
                    state = 'fragile'
                else:
                    state = 'absent'
                grid[fi, pi * len(radii) + ri] = value_map[state]
    fig, ax = plt.subplots(figsize=(12.4, 4.2))
    im = ax.imshow(grid, aspect='auto', cmap='viridis', vmin=0.0, vmax=2.0)
    xticks = []
    for phase in phases:
        for radius in radii:
            xticks.append(f'{phase:.3f}\n{radius:+d}')
    ax.set_xticks(range(len(xticks)), xticks, fontsize=7)
    ax.set_yticks(range(len(families)), ['clustered', 'corridor', 'triad'])
    ax.set_title('Stage 23.9 topology stability map (family x phase x radius)')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_midpoint_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    clustered = [row for row in rows if row['family_id'] == 'clustered_dk_baseline' and row['kernel_radius_offset'] == 0 and row['resolution'] == 12]
    clustered.sort(key=lambda row: float(row['phase_detuning_fraction_of_pi']))
    phases = [float(row['phase_detuning_fraction_of_pi']) for row in clustered]
    flow = [float(row['flow_concentration_index']) for row in clustered]
    topo_map = {'braid_like_exchange': 0, 'transfer_smeared': 1, 'localized_capture': 2, 'unresolved_mixed': 3}
    labels = [topo_map[str(row['topology_class'])] for row in clustered]
    fig, axes = plt.subplots(2, 1, figsize=(6.8, 6.0), sharex=True)
    axes[0].plot(phases, flow, marker='o', color='tab:blue')
    axes[0].set_ylabel('flow concentration')
    axes[0].set_title('Stage 23.9 midpoint corridor characterization')
    axes[0].grid(alpha=0.25)
    axes[1].plot(phases, labels, marker='o', color='tab:orange')
    axes[1].set_yticks([0, 1, 2, 3], ['braid', 'smeared', 'capture', 'mixed'])
    axes[1].set_xlabel('phase detuning (fraction of pi)')
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_transport_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    families = ['clustered_dk_baseline', 'corridor_transport_family', 'triad_competition_family']
    counts = [sum(int(row['metastable_transport_flag']) for row in rows if row['family_id'] == family) for family in families]
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.bar(range(len(families)), counts, color='tab:red')
    ax.set_xticks(range(len(families)), ['clustered', 'corridor', 'triad'])
    ax.set_ylabel('metastable transport hits')
    ax.set_title('Stage 23.9 transport-window confirmation')
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_refinement_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    families = ['clustered_dk_baseline', 'corridor_transport_family', 'triad_competition_family']
    vals = [float(np.mean([float(row['refinement_survival_ratio']) for row in rows if row['family_id'] == family])) for family in families]
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.bar(range(len(families)), vals, color='tab:green')
    ax.set_xticks(range(len(families)), ['clustered', 'corridor', 'triad'])
    ax.set_ylabel('mean refinement survival ratio')
    ax.set_title('Stage 23.9 refinement consistency')
    ax.set_ylim(0.0, 1.05)
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_irreversibility_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    families = ['clustered_dk_baseline', 'corridor_transport_family', 'triad_competition_family']
    vals = [float(np.max([float(row['structural_irreversibility_score']) for row in rows if row['family_id'] == family])) for family in families]
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.bar(range(len(families)), vals, color='tab:purple')
    ax.set_xticks(range(len(families)), ['clustered', 'corridor', 'triad'])
    ax.set_ylabel('max irreversibility score')
    ax.set_title('Stage 23.9 limited irreversibility mapping')
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    json_path: Path,
    csv_path: Path,
    rows: list[dict[str, Any]],
    stamped_plots: list[str],
    outcome: str,
) -> None:
    topo_counts = Counter(str(row['topology_class']) for row in rows)
    clustered = [row for row in rows if row['family_id'] == 'clustered_dk_baseline']
    corridor = [row for row in rows if row['family_id'] == 'corridor_transport_family']
    triad = [row for row in rows if row['family_id'] == 'triad_competition_family']
    midpoint_hits = sum(int(row['midpoint_band_membership_flag']) for row in clustered)
    transport_hits = sum(int(row['metastable_transport_flag']) for row in rows)
    lines = [
        '# Stage 23.9 DK Collision Phenomenology Consolidation v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Outcome: `{outcome}`',
        '',
        '## Topology counts',
    ]
    for label, count in sorted(topo_counts.items()):
        lines.append(f'- `{label}`: {count}')
    lines.extend([
        '',
        '## Consolidation read',
        f"- clustered braid hits: {sum(int(row['braid_like_exchange']) for row in clustered)} / {len(clustered)}",
        f"- midpoint-band clustered smear hits: {midpoint_hits}",
        f"- transport-window hits: {transport_hits}",
        f"- corridor localized-capture hits: {sum(int(row['localized_capture']) for row in corridor)} / {len(corridor)}",
        f"- triad braid hits: {sum(int(row['braid_like_exchange']) for row in triad)} / {len(triad)}",
        f"- max family irreversibility scores: clustered={max(float(row['structural_irreversibility_score']) for row in clustered):.4f}, corridor={max(float(row['structural_irreversibility_score']) for row in corridor):.4f}, triad={max(float(row['structural_irreversibility_score']) for row in triad):.4f}",
        '',
        '## Clean conclusion',
    ])
    if outcome == 'phenomenology_consolidated':
        lines.append('Phase III now supports a stable phenomenology layer: the clustered family retains a reproducible braid-exchange window, the smeared midpoint corridor remains finite, no metastable transport family opens, and limited irreversibility pockets remain bounded and seed-dependent rather than persistence-level.')
    elif outcome == 'partial_structural_instability':
        lines.append('Phase III retains the main clustered braid phenomenology and the transport null, but at least one family or axis remains too unstable to freeze the whole lattice as a fully consolidated atlas layer.')
    else:
        lines.append('The frozen operator branch does not stabilize into a coherent common lattice across the tested families, so the present pass is family-inconsistent rather than fully consolidated.')
    lines.extend([
        '',
        '## Plots',
    ])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')
    NOTE_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    common = runsheet['common_fields']

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for family in runsheet['families']:
        for phase_fraction in common['phase_detunings_fraction_of_pi']:
            for radius_offset in common['kernel_radius_offsets_in_steps']:
                for resolution in (int(common['base_resolution']), int(common['refined_resolution'])):
                    row, detail, run_plots = simulate_run(
                        family=family,
                        phase_fraction=float(phase_fraction),
                        radius_offset=int(radius_offset),
                        resolution=int(resolution),
                        common=common,
                    )
                    rows.append(row)
                    details.append(detail)
                    plot_paths.extend(run_plots)

    apply_refinement_survival(rows)
    rows.sort(key=lambda row: row['run_id'])
    details.sort(key=lambda detail: detail['metrics']['run_id'])

    stability_map = WORK_PLOT_DIR / 'stage23_9_topology_stability_map.png'
    midpoint_panel = WORK_PLOT_DIR / 'stage23_9_midpoint_corridor_panel.png'
    transport_panel = WORK_PLOT_DIR / 'stage23_9_transport_window_panel.png'
    refinement_panel = WORK_PLOT_DIR / 'stage23_9_refinement_consistency_panel.png'
    irreversibility_panel = WORK_PLOT_DIR / 'stage23_9_irreversibility_panel.png'
    plot_topology_stability_map(stability_map, rows)
    plot_midpoint_panel(midpoint_panel, rows)
    plot_transport_panel(transport_panel, rows)
    plot_refinement_panel(refinement_panel, rows)
    plot_irreversibility_panel(irreversibility_panel, rows)
    plot_paths.extend([stability_map, midpoint_panel, transport_panel, refinement_panel, irreversibility_panel])

    outcome = classify_outcome(rows)
    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': {
            'outcome': outcome,
            'topology_counts': dict(Counter(str(row['topology_class']) for row in rows)),
            'transport_window_hits': int(sum(int(row['metastable_transport_flag']) for row in rows)),
            'clustered_midpoint_band_hits': int(sum(int(row['midpoint_band_membership_flag']) for row in rows if row['family_id'] == 'clustered_dk_baseline')),
        },
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage23_9_dk_collision_phenomenology_consolidation',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, outcome)


if __name__ == '__main__':
    main()
