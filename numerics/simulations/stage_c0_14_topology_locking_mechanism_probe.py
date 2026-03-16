#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
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
    connected_components,
    edge_grid,
    evolve,
    grade_exchange_signal,
    grade_weights,
    pair_phase_difference,
    recurrence_indicator,
    packet_state_with_skew,
)
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions
from stage_c0_4_controlled_motif_injection import choose_seed_patch, triangle_specs
from stage_c0_11_harmonic_address_resonance_ladder import combined_initial_weights, operator_adjacency_lists
from stage_c0_12_harmonic_detuning_continuum_scan import (
    address_activity_statistics,
    address_weighted_operator,
    build_address_labels,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_14_topology_locking_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_14_Topology_Locking_Mechanism_Probe_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'run_order',
    'family_id',
    'family_label',
    'family_mode',
    'lock_rule',
    'lock_label',
    'delta',
    'topology_class',
    'locking_label',
    'braid_like_exchange',
    'transfer_smeared',
    'locked_braid',
    'quasi_bounded_exchange',
    'unresolved',
    'topology_survival_time',
    'baseline_topology_survival_time',
    'survival_time_gain',
    'recurrence_indicator',
    'local_path_reuse_score',
    'harmonic_latch_duty_cycle',
    'motif_anchor_score',
    'flow_concentration_index',
    'exchange_coherence',
    'transport_span',
    'bounded_dwell_proxy',
    'baseline_topology_class',
    'baseline_flow_concentration_index',
    'baseline_exchange_coherence',
    'baseline_transport_span',
    'baseline_bounded_dwell_proxy',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.14 topology-locking mechanism probe.')
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


def make_rng(name: str) -> np.random.Generator:
    seed = 5003 + sum((idx + 1) * ord(ch) for idx, ch in enumerate(name))
    return np.random.default_rng(seed)


def build_family_labels_with_mode(
    packet_states: list[np.ndarray],
    block_sizes: tuple[int, ...],
    delta: float,
    family_mode: str,
    support_floor_fraction: float,
    dominant_ratio_threshold: float,
    modulus: float,
    rng: np.random.Generator,
) -> np.ndarray:
    left_address = 0.0
    right_address = 2.0 * float(delta)
    if family_mode == 'randomized':
        labels = rng.uniform(0.0, float(modulus), size=sum(block_sizes))
        return np.asarray(labels, dtype=float)
    return build_address_labels(
        packet_states=packet_states,
        block_sizes=block_sizes,
        left_address=left_address,
        right_address=right_address,
        support_floor_fraction=support_floor_fraction,
        dominant_ratio_threshold=dominant_ratio_threshold,
    )


def node_profile(states: list[np.ndarray]) -> np.ndarray:
    weights = np.mean(np.asarray([np.abs(state) ** 2 for state in states], dtype=float), axis=0)
    total = float(np.sum(weights))
    if total <= 1.0e-12:
        return np.zeros_like(weights)
    return weights / total


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 1.0e-12:
        return 0.0
    return float(np.clip(np.dot(a, b) / denom, 0.0, 1.0))


def apply_node_bias(operator: sp.csr_matrix, node_bias: np.ndarray, amplitude: float) -> sp.csr_matrix:
    if amplitude <= 0.0 or not np.any(node_bias > 1.0e-12):
        return operator
    scaled = np.asarray(node_bias, dtype=float)
    peak = float(np.max(scaled))
    if peak <= 1.0e-12:
        return operator
    scaled = scaled / peak
    coo = operator.tocoo()
    mult = 1.0 + float(amplitude) * 0.5 * (scaled[coo.row] + scaled[coo.col])
    return sp.csr_matrix((coo.data * mult, (coo.row, coo.col)), shape=operator.shape).tocsr()


def apply_harmonic_latch_bias(
    operator: sp.csr_matrix,
    labels: np.ndarray,
    modulus: float,
    latch_level: float,
    amplitude: float,
) -> sp.csr_matrix:
    if latch_level <= 1.0e-12 or amplitude <= 0.0:
        return operator
    coo = operator.tocoo()
    diff = np.mod(np.abs(labels[coo.row] - labels[coo.col]), modulus)
    distance = np.minimum(diff, modulus - diff)
    compat = np.clip(1.0 - 0.5 * distance, 0.0, 1.0)
    mult = 1.0 + float(amplitude) * float(latch_level) * compat
    return sp.csr_matrix((coo.data * mult, (coo.row, coo.col)), shape=operator.shape).tocsr()


def identify_motif_nodes(
    base_operator: sp.csr_matrix,
    packet_states: list[np.ndarray],
    params: dict[str, Any],
) -> list[int]:
    adjacency_lists = operator_adjacency_lists(base_operator)
    avg_weights = combined_initial_weights(packet_states)
    specs = triangle_specs(adjacency_lists, avg_weights, params)
    if specs:
        return [int(value) for value in specs[0]['motif_nodes']]
    _seed, _patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=3, params=params)
    return [int(value) for value in candidates[:3]]


def measure_window(
    states: list[np.ndarray],
    packet_states: list[np.ndarray],
    positions: np.ndarray,
    edge_midpoints: np.ndarray,
    block_sizes: tuple[int, ...],
    resolution: int,
    sigma: float,
    labels: np.ndarray,
    operator: sp.csr_matrix,
    modulus: float,
    eta_match: float,
    eta_mismatch: float,
) -> dict[str, Any]:
    n0, n1, _ = block_sizes
    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(2)]
    previous: list[np.ndarray] | None = None
    braid_votes: list[float] = []
    edge_grids: list[np.ndarray] = []
    phase_diffs: list[float] = []
    raw_peak_counts: list[int] = []
    for idx, state in enumerate(states):
        grid = field_grid_2d(positions, state, resolution)
        current, raw_count = ordered_peaks(grid, 2, None if idx else None, previous)
        raw_peak_counts.append(raw_count)
        for track, point in zip(peak_tracks, current):
            track.append(point.copy())
        if len(current) == 2:
            braid_votes.append(float(current[0][0] > current[1][0]))
        previous = [point.copy() for point in current]
        edge_values = np.abs(state[n0:n0 + n1])
        edge_grids.append(edge_grid(edge_midpoints, edge_values, resolution))
        phase_diffs.append(pair_phase_difference(state, packet_states))

    center_histories = [np.asarray(track, dtype=float) for track in peak_tracks if track]
    if len(center_histories) == 2:
        mean_pair_distances, min_pair_distances = mean_pair_separation_series(center_histories)
    else:
        mean_pair_distances = np.asarray([], dtype=float)
        min_pair_distances = np.asarray([], dtype=float)
    close_threshold = 2.0 * sigma
    avg_edge_grid = np.mean(np.asarray(edge_grids, dtype=float), axis=0) if edge_grids else np.zeros((resolution, resolution), dtype=float)
    total = float(np.sum(avg_edge_grid))
    flat = np.sort(avg_edge_grid.ravel())[::-1]
    top_k = max(1, int(math.ceil(0.1 * flat.size)))
    concentration = float(np.sum(flat[:top_k]) / max(total, 1.0e-12))
    grade_hist = np.asarray([grade_weights(state, block_sizes) for state in states], dtype=float)
    transfer_signal, coherence = grade_exchange_signal(grade_hist)
    _same_activity, sorting, selectivity = address_activity_statistics(states, operator, labels, modulus, eta_match, eta_mismatch)
    braid_flag = abs(float(np.mean(braid_votes)) - 0.5) > 0.25 if braid_votes else False
    topology = classify_selection_topology(concentration, coherence, braid_flag, selectivity)
    recurrence = recurrence_indicator(mean_pair_distances, close_threshold) if mean_pair_distances.size else 0.0
    survival = float(np.sum(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    final_mean = float(mean_pair_distances[-1]) if mean_pair_distances.size else 0.0
    min_sep = float(np.min(min_pair_distances)) if min_pair_distances.size else 0.0
    transport_span = float(max(0.0, final_mean - min_sep))
    bounded_dwell = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    threshold = max(0.15 * float(np.max(avg_edge_grid)), float(np.mean(avg_edge_grid) + np.std(avg_edge_grid)))
    channels = connected_components(avg_edge_grid >= threshold)
    return {
        'avg_edge_grid': avg_edge_grid,
        'flow_concentration_index': concentration,
        'exchange_coherence': coherence,
        'address_selectivity_index': selectivity,
        'address_sorting_score': sorting,
        'braid_flag': braid_flag,
        'topology_class': topology,
        'recurrence_indicator': recurrence,
        'survival_count': survival,
        'transport_span': transport_span,
        'bounded_dwell_proxy': bounded_dwell,
        'local_channel_count': int(channels),
        'phase_alignment_metric': phase_alignment_metric_from_list(phase_diffs),
        'raw_peak_count': int(raw_peak_counts[-1]) if raw_peak_counts else 0,
        'close_threshold': close_threshold,
        'mean_pair_distances': mean_pair_distances,
        'transfer_signal': transfer_signal,
    }


def classify_selection_topology(concentration: float, coherence: float, braid_flag: bool, selectivity: float) -> str:
    if braid_flag and concentration >= 0.89 and coherence >= 0.44 and selectivity >= 0.88:
        return 'braid_like_exchange'
    if concentration >= 0.84 and coherence >= 0.34:
        return 'transfer_smeared'
    return 'unresolved'


def phase_alignment_metric_from_list(phase_diffs: list[float]) -> float:
    if not phase_diffs:
        return 0.0
    mean_abs = float(np.mean(np.abs(np.asarray(phase_diffs, dtype=float))))
    return float(np.clip(1.0 - mean_abs / math.pi, 0.0, 1.0))


def distributed_segments(total_steps: int, segment_count: int) -> list[int]:
    segment_count = max(1, int(segment_count))
    base = total_steps // segment_count
    remainder = total_steps % segment_count
    out = []
    for idx in range(segment_count):
        out.append(base + (1 if idx < remainder else 0))
    return [max(1, value) for value in out]


def topology_code(label: str) -> float:
    return {
        'braid_like_exchange': 0.0,
        'transfer_smeared': 1.0,
        'locked_braid': 2.0,
        'quasi_bounded_exchange': 3.0,
        'unresolved': 4.0,
    }[label]


def locking_rule_code(label: str) -> float:
    return {
        'no_locking_effect': 0.0,
        'topology_stabilized_no_capture': 1.0,
        'quasi_bounded_exchange': 2.0,
        'motif_localized_lock': 3.0,
        'address_latched_exchange': 4.0,
    }[label]


def final_topology_class(
    selection_topology: str,
    survival_gain: float,
    recurrence: float,
    path_reuse: float,
    latch_duty: float,
    motif_anchor: float,
    transport_span: float,
    baseline_transport_span: float,
    bounded_dwell: float,
    baseline_bounded_dwell: float,
) -> str:
    if (
        selection_topology == 'braid_like_exchange'
        and survival_gain >= 0.10
        and bounded_dwell >= baseline_bounded_dwell + 0.08
        and transport_span <= baseline_transport_span + 0.02
        and max(recurrence, path_reuse, latch_duty, motif_anchor) >= 0.24
    ):
        return 'locked_braid'
    if (
        survival_gain >= 0.06
        and bounded_dwell >= baseline_bounded_dwell + 0.04
        and transport_span <= baseline_transport_span + 0.05
        and max(path_reuse, latch_duty, motif_anchor) >= 0.18
    ):
        return 'quasi_bounded_exchange'
    if selection_topology == 'braid_like_exchange':
        return 'braid_like_exchange'
    if selection_topology == 'transfer_smeared':
        return 'transfer_smeared'
    return 'unresolved'


def derived_locking_label(
    lock_rule: str,
    final_topology: str,
    survival_gain: float,
    path_reuse: float,
    latch_duty: float,
    motif_anchor: float,
) -> str:
    if lock_rule == 'harmonic_phase_latch' and final_topology in {'locked_braid', 'quasi_bounded_exchange'} and latch_duty >= 0.35:
        return 'address_latched_exchange'
    if lock_rule == 'motif_occupancy_latch' and final_topology in {'locked_braid', 'quasi_bounded_exchange'} and motif_anchor >= 0.12:
        return 'motif_localized_lock'
    if final_topology == 'quasi_bounded_exchange':
        return 'quasi_bounded_exchange'
    if final_topology == 'locked_braid' or (
        survival_gain >= 0.08 and final_topology in {'braid_like_exchange', 'quasi_bounded_exchange'}
    ):
        return 'topology_stabilized_no_capture'
    return 'no_locking_effect'


def summarize_measurement(measurement: dict[str, Any]) -> dict[str, Any]:
    return {
        'topology_class': str(measurement['topology_class']),
        'flow_concentration_index': float(measurement['flow_concentration_index']),
        'exchange_coherence': float(measurement['exchange_coherence']),
        'address_selectivity_index': float(measurement['address_selectivity_index']),
        'address_sorting_score': float(measurement['address_sorting_score']),
        'recurrence_indicator': float(measurement['recurrence_indicator']),
        'transport_span': float(measurement['transport_span']),
        'bounded_dwell_proxy': float(measurement['bounded_dwell_proxy']),
        'local_channel_count': int(measurement['local_channel_count']),
        'phase_alignment_metric': float(measurement['phase_alignment_metric']),
        'raw_peak_count': int(measurement['raw_peak_count']),
        'survival_count': float(measurement['survival_count']),
    }


def topology_trace_duration(trace: list[str], segment_lengths: list[int], dt: float) -> float:
    total = 0.0
    for label, steps in zip(trace, segment_lengths):
        if label == 'braid_like_exchange':
            total += float(steps) * float(dt)
    return total


def simulate_family(
    run_id: str,
    family_mode: str,
    delta: float,
    lock_rule: str,
    common: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
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
    segment_count = int(common['segment_count'])

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes

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

    rng = make_rng(run_id)
    labels = build_family_labels_with_mode(
        packet_states=packet_states,
        block_sizes=block_sizes,
        delta=delta,
        family_mode=family_mode,
        support_floor_fraction=support_floor_fraction,
        dominant_ratio_threshold=dominant_ratio_threshold,
        modulus=modulus,
        rng=rng,
    )
    base_operator = address_weighted_operator(complex_data.dirac_kahler, labels, modulus, eta_match, eta_mismatch)
    motif_nodes = identify_motif_nodes(base_operator, packet_states, common)

    dt = first_order_dt(complex_data.dirac_kahler, dt_scale)
    steps = max(2, int(math.ceil(t_final / dt)))
    segment_lengths = distributed_segments(steps, segment_count)
    segment_times = np.cumsum(np.asarray(segment_lengths, dtype=float)) * dt

    def run_sequence(active_lock_rule: str) -> dict[str, Any]:
        psi = np.sum(packet_states, axis=0)
        all_states = [psi.copy()]
        path_latch = np.zeros_like(psi, dtype=float)
        harmonic_latch = 0.0
        motif_latch = 0.0
        prev_profile = None
        path_reuse_scores: list[float] = []
        latch_scores: list[float] = []
        motif_scores: list[float] = []
        flow_trace: list[float] = []
        topology_trace: list[str] = []
        topology_trace_values: list[float] = []
        for seg_idx, seg_steps in enumerate(segment_lengths):
            operator = base_operator
            if active_lock_rule == 'exchange_path_reinforcement':
                operator = apply_node_bias(operator, path_latch, float(common['path_lock_amplitude']))
            elif active_lock_rule == 'harmonic_phase_latch':
                operator = apply_harmonic_latch_bias(operator, labels, modulus, harmonic_latch, float(common['harmonic_latch_amplitude']))
            elif active_lock_rule == 'motif_occupancy_latch':
                motif_bias = np.zeros_like(psi, dtype=float)
                motif_bias[motif_nodes] = motif_latch
                operator = apply_node_bias(operator, motif_bias, float(common['motif_latch_amplitude']))

            segment_states = evolve(operator, block_sizes, psi, dt, int(seg_steps), beta)
            psi = segment_states[-1]
            all_states.extend(segment_states[1:])

            window = measure_window(
                states=segment_states,
                packet_states=packet_states,
                positions=positions,
                edge_midpoints=complex_data.edge_midpoints,
                block_sizes=block_sizes,
                resolution=resolution,
                sigma=sigma,
                labels=labels,
                operator=operator,
                modulus=modulus,
                eta_match=eta_match,
                eta_mismatch=eta_mismatch,
            )
            profile = node_profile(segment_states)
            reuse = cosine_similarity(profile, prev_profile) if prev_profile is not None else 0.0
            prev_profile = profile
            path_reuse_scores.append(reuse)

            flow_trace.append(float(window['flow_concentration_index']))
            topology_trace.append(str(window['topology_class']))
            topology_trace_values.append(topology_code(str(window['topology_class']).replace('address_protected_braid', 'braid_like_exchange') if 'address_protected_braid' in str(window['topology_class']) else str(window['topology_class'])))

            if active_lock_rule == 'exchange_path_reinforcement':
                trigger = 1.0 if str(window['topology_class']) == 'braid_like_exchange' else 0.5 if float(window['flow_concentration_index']) >= 0.86 else 0.0
                path_latch = float(common['path_lock_decay']) * path_latch + trigger * profile
                peak = float(np.max(path_latch))
                if peak > 1.0e-12:
                    path_latch = path_latch / peak
                latch_scores.append(float(np.mean(path_latch[profile > np.mean(profile)])) if np.any(profile > np.mean(profile)) else 0.0)
                motif_scores.append(0.0)
            elif active_lock_rule == 'harmonic_phase_latch':
                trigger = 1.0 if float(window['address_selectivity_index']) >= 0.88 and float(window['exchange_coherence']) >= 0.42 else 0.0
                harmonic_latch = min(1.0, float(common['harmonic_latch_decay']) * harmonic_latch + 0.55 * trigger)
                latch_scores.append(harmonic_latch)
                motif_scores.append(0.0)
            elif active_lock_rule == 'motif_occupancy_latch':
                occupancy = float(np.sum(profile[motif_nodes])) if motif_nodes else 0.0
                trigger = occupancy if str(window['topology_class']) == 'braid_like_exchange' else 0.5 * occupancy
                motif_latch = min(1.0, float(common['motif_latch_decay']) * motif_latch + trigger)
                latch_scores.append(motif_latch)
                motif_scores.append(occupancy)
            else:
                latch_scores.append(0.0)
                motif_scores.append(float(np.sum(profile[motif_nodes])) if motif_nodes else 0.0)

        full = measure_window(
            states=all_states,
            packet_states=packet_states,
            positions=positions,
            edge_midpoints=complex_data.edge_midpoints,
            block_sizes=block_sizes,
            resolution=resolution,
            sigma=sigma,
            labels=labels,
            operator=base_operator,
            modulus=modulus,
            eta_match=eta_match,
            eta_mismatch=eta_mismatch,
        )
        return {
            'full': full,
            'times': [idx * dt for idx in range(len(all_states))],
            'segment_times': segment_times.tolist(),
            'all_states': all_states,
            'topology_trace': topology_trace,
            'topology_trace_values': topology_trace_values,
            'flow_trace': flow_trace,
            'path_reuse_scores': path_reuse_scores,
            'latch_scores': latch_scores,
            'motif_scores': motif_scores,
            'motif_nodes': motif_nodes,
        }

    baseline = run_sequence('none')
    locked = run_sequence(lock_rule)

    baseline_full = baseline['full']
    locked_full = locked['full']
    topology_survival_time = min(topology_trace_duration(locked['topology_trace'], segment_lengths, dt), t_final)
    baseline_survival = min(topology_trace_duration(baseline['topology_trace'], segment_lengths, dt), t_final)
    survival_gain = topology_survival_time - baseline_survival
    path_reuse = float(np.mean(locked['path_reuse_scores'])) if locked['path_reuse_scores'] else 0.0
    latch_duty = float(np.mean(np.asarray(locked['latch_scores']) > 0.25)) if locked['latch_scores'] else 0.0
    motif_anchor = float(np.mean(locked['motif_scores'])) if locked['motif_scores'] else 0.0
    final_topology = final_topology_class(
        selection_topology=str(locked_full['topology_class']),
        survival_gain=survival_gain,
        recurrence=float(locked_full['recurrence_indicator']),
        path_reuse=path_reuse,
        latch_duty=latch_duty,
        motif_anchor=motif_anchor,
        transport_span=float(locked_full['transport_span']),
        baseline_transport_span=float(baseline_full['transport_span']),
        bounded_dwell=float(locked_full['bounded_dwell_proxy']),
        baseline_bounded_dwell=float(baseline_full['bounded_dwell_proxy']),
    )
    locking_label = derived_locking_label(lock_rule, final_topology, survival_gain, path_reuse, latch_duty, motif_anchor)

    row = {
        'run_id': run_id,
        'geometry_id': 'tight_clustered_pair',
        'phase_id': f'lock_{lock_rule}',
        'resolution': resolution,
        'graph_type': common['graph_type'],
        'packet_count': 2,
        'initial_peak_count': int(locked_full['raw_peak_count']),
        'collision_label': final_topology,
        'persistence_label': locking_label,
        'composite_lifetime': topology_survival_time,
        'binding_persistence': float(locked_full['bounded_dwell_proxy']),
        'corridor_dwell': 0.0,
        'coarse_basin_persistence': 0.0,
        'minimum_separation': float(np.min(locked_full['mean_pair_distances'])) if locked_full['mean_pair_distances'].size else 0.0,
        'final_mean_separation': float(locked_full['mean_pair_distances'][-1]) if locked_full['mean_pair_distances'].size else 0.0,
        'post_collision_separation_trend': float(locked_full['transport_span']),
        'encounter_dwell_time': topology_survival_time,
        'deflection_angle_proxy': 0.0,
        'reflection_fraction': 0.0,
        'grade_transfer_amplitude': float(np.max(locked_full['transfer_signal'])) if locked_full['transfer_signal'].size else 0.0,
        'omega0_weight_initial': 0.0,
        'omega1_weight_initial': 0.0,
        'omega2_weight_initial': 0.0,
        'omega0_weight_final': 0.0,
        'omega1_weight_final': 0.0,
        'omega2_weight_final': 0.0,
        'localized_circulation_proxy': 0.0,
        'new_collision_class': int(final_topology in {'locked_braid', 'quasi_bounded_exchange'}),
        'gate_met': int(locking_label != 'no_locking_effect'),
        'promoted_followup': 0,
        'run_order': 0,
        'family_id': '',
        'family_label': '',
        'family_mode': family_mode,
        'lock_rule': lock_rule,
        'lock_label': '',
        'delta': float(delta),
        'topology_class': final_topology,
        'locking_label': locking_label,
        'braid_like_exchange': int(final_topology == 'braid_like_exchange'),
        'transfer_smeared': int(final_topology == 'transfer_smeared'),
        'locked_braid': int(final_topology == 'locked_braid'),
        'quasi_bounded_exchange': int(final_topology == 'quasi_bounded_exchange'),
        'unresolved': int(final_topology == 'unresolved'),
        'topology_survival_time': topology_survival_time,
        'baseline_topology_survival_time': baseline_survival,
        'survival_time_gain': survival_gain,
        'recurrence_indicator': float(locked_full['recurrence_indicator']),
        'local_path_reuse_score': path_reuse,
        'harmonic_latch_duty_cycle': latch_duty,
        'motif_anchor_score': motif_anchor,
        'flow_concentration_index': float(locked_full['flow_concentration_index']),
        'exchange_coherence': float(locked_full['exchange_coherence']),
        'transport_span': float(locked_full['transport_span']),
        'bounded_dwell_proxy': float(locked_full['bounded_dwell_proxy']),
        'baseline_topology_class': str(baseline_full['topology_class']),
        'baseline_flow_concentration_index': float(baseline_full['flow_concentration_index']),
        'baseline_exchange_coherence': float(baseline_full['exchange_coherence']),
        'baseline_transport_span': float(baseline_full['transport_span']),
        'baseline_bounded_dwell_proxy': float(baseline_full['bounded_dwell_proxy']),
        'notes': '',
    }

    detail = {
        'baseline': {
            'summary': summarize_measurement(baseline_full),
            'segment_times': baseline['segment_times'],
            'topology_trace': baseline['topology_trace'],
            'flow_trace': baseline['flow_trace'],
            'path_reuse_scores': baseline['path_reuse_scores'],
            'latch_scores': baseline['latch_scores'],
            'motif_scores': baseline['motif_scores'],
        },
        'locked': {
            'summary': summarize_measurement(locked_full),
            'segment_times': locked['segment_times'],
            'topology_trace': locked['topology_trace'],
            'flow_trace': locked['flow_trace'],
            'path_reuse_scores': locked['path_reuse_scores'],
            'latch_scores': locked['latch_scores'],
            'motif_scores': locked['motif_scores'],
            'motif_nodes': [int(value) for value in locked['motif_nodes']],
        },
        'metrics': row,
    }

    plot_paths: list[Path] = []
    topology_path = WORK_PLOT_DIR / f"stage_c0_14_{run_id}_topology_trace.png"
    latch_path = WORK_PLOT_DIR / f"stage_c0_14_{run_id}_lock_score_trace.png"
    flow_path = WORK_PLOT_DIR / f"stage_c0_14_{run_id}_flow_concentration_trace.png"
    plot_trace(topology_path, run_id, locked['segment_times'], np.asarray(locked['topology_trace_values'], dtype=float), 'topology code', 'Stage C0.14 topology trace')
    plot_trace(latch_path, run_id, locked['segment_times'], np.asarray(locked['latch_scores'], dtype=float), 'lock score', 'Stage C0.14 latch / lock score trace')
    plot_trace(flow_path, run_id, locked['segment_times'], np.asarray(locked['flow_trace'], dtype=float), 'flow concentration', 'Stage C0.14 flow concentration trace')
    plot_paths.extend([topology_path, latch_path, flow_path])
    return row, detail, plot_paths


def plot_rule_family_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    families = [str(row['family_label']) for row in sorted(rows, key=lambda item: int(item['run_order'])) if str(row['lock_rule']) == 'exchange_path_reinforcement']
    rules = ['exchange_path_reinforcement', 'harmonic_phase_latch', 'motif_occupancy_latch']
    rule_labels = ['Lock A', 'Lock B', 'Lock C']
    matrix = np.zeros((len(rules), len(families)), dtype=float)
    texts: dict[tuple[int, int], str] = {}
    for i, rule in enumerate(rules):
        for j, family in enumerate(families):
            row = next(item for item in rows if str(item['lock_rule']) == rule and str(item['family_label']) == family)
            matrix[i, j] = locking_rule_code(str(row['locking_label']))
            texts[(i, j)] = str(row['topology_class']).replace('_', ' ')
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    im = ax.imshow(matrix, aspect='auto', cmap='viridis', vmin=0.0, vmax=4.0)
    ax.set_xticks(range(len(families)), families)
    ax.set_yticks(range(len(rule_labels)), rule_labels)
    ax.set_title('Stage C0.14 locking-rule x family class matrix')
    for (i, j), text in texts.items():
        ax.text(j, i, text, ha='center', va='center', color='white', fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_survival_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [f"{row['family_label']}\n{row['lock_label']}" for row in rows]
    baseline = np.asarray([float(row['baseline_topology_survival_time']) for row in rows], dtype=float)
    locked = np.asarray([float(row['topology_survival_time']) for row in rows], dtype=float)
    x = np.arange(len(rows))
    width = 0.38
    fig, ax = plt.subplots(figsize=(11.2, 5.0))
    ax.bar(x - width / 2.0, baseline, width=width, label='baseline')
    ax.bar(x + width / 2.0, locked, width=width, label='locked rule')
    ax.set_xticks(x, labels, rotation=35, ha='right')
    ax.set_ylabel('topology survival time')
    ax.set_title('Stage C0.14 topology survival comparison')
    ax.legend()
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_latch_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [f"{row['family_label']}\n{row['lock_label']}" for row in rows]
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(11.2, 5.0))
    ax.plot(x, [float(row['local_path_reuse_score']) for row in rows], marker='o', label='path reuse')
    ax.plot(x, [float(row['harmonic_latch_duty_cycle']) for row in rows], marker='s', label='harmonic latch duty')
    ax.plot(x, [float(row['motif_anchor_score']) for row in rows], marker='^', label='motif anchor')
    ax.set_xticks(x, labels, rotation=35, ha='right')
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel('score')
    ax.set_title('Stage C0.14 latch / anchor summary')
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def global_verdict(rows: list[dict[str, Any]]) -> str:
    positive = [
        row for row in rows
        if str(row['locking_label']) in {'address_latched_exchange', 'motif_localized_lock', 'quasi_bounded_exchange'}
        or str(row['topology_class']) in {'locked_braid', 'quasi_bounded_exchange'}
    ]
    stabilized = [row for row in rows if str(row['locking_label']) == 'topology_stabilized_no_capture']
    if positive:
        return 'at least one bounded combinatorial latch converts topology selection into a longer-lived exchange structure'
    if stabilized:
        return 'locking rules extend topology survival in places, but no new locking class appears'
    return 'topology selection persists, but the tested combinatorial locking rules do not create locking'


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], verdict: str) -> None:
    lines = [
        '# Stage C0.14 Topology Locking Mechanism Probe v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Global read: `{verdict}`',
        '',
        '## Locking-rule x family summary',
    ]
    for row in rows:
        lines.append(
            f"- `{row['run_id']}`: family=`{row['family_label']}`, lock=`{row['lock_label']}`, topology=`{row['topology_class']}`, locking_label=`{row['locking_label']}`, survival_gain=`{row['survival_time_gain']:.4f}`, recurrence=`{row['recurrence_indicator']:.4f}`"
        )
    lines.extend([
        '',
        '## Score summary',
    ])
    for row in rows:
        lines.append(
            f"- `{row['run_id']}`: path_reuse=`{row['local_path_reuse_score']:.4f}`, latch_duty=`{row['harmonic_latch_duty_cycle']:.4f}`, motif_anchor=`{row['motif_anchor_score']:.4f}`, bounded_dwell=`{row['bounded_dwell_proxy']:.4f}`"
        )
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
    runs = selected_runs(runsheet['runs'], args.run_ids)
    common = runsheet['common_fields']

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        row, detail, run_plots = simulate_family(
            run_id=str(run['run_id']),
            family_mode=str(run['family_mode']),
            delta=float(run['delta']),
            lock_rule=str(run['lock_rule']),
            common=common,
        )
        row['run_order'] = int(run['run_order'])
        row['family_id'] = str(run['family_id'])
        row['family_label'] = str(run['family_label'])
        row['family_mode'] = str(run['family_mode'])
        row['lock_rule'] = str(run['lock_rule'])
        row['lock_label'] = str(run['lock_label'])
        row['notes'] = str(run['notes'])
        detail['metrics'] = row
        detail['run'] = run
        rows.append(row)
        details.append(detail)
        plot_paths.extend(run_plots)

    rows.sort(key=lambda item: int(item['run_order']))
    details.sort(key=lambda item: int(item['metrics']['run_order']))

    matrix_path = WORK_PLOT_DIR / 'stage_c0_14_locking_rule_family_matrix.png'
    survival_path = WORK_PLOT_DIR / 'stage_c0_14_survival_comparison_panel.png'
    latch_path = WORK_PLOT_DIR / 'stage_c0_14_latch_anchor_summary.png'
    plot_rule_family_matrix(matrix_path, rows)
    plot_survival_panel(survival_path, rows)
    plot_latch_summary(latch_path, rows)
    plot_paths.extend([matrix_path, survival_path, latch_path])

    verdict = global_verdict(rows)
    summary = {
        'topology_counts': dict(Counter(str(row['topology_class']) for row in rows)),
        'locking_counts': dict(Counter(str(row['locking_label']) for row in rows)),
        'mean_survival_gain_by_rule': {
            rule: float(np.mean([float(row['survival_time_gain']) for row in rows if str(row['lock_rule']) == rule]))
            for rule in sorted({str(row['lock_rule']) for row in rows})
        },
        'verdict': verdict,
    }
    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': summary,
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage_c0_14_topology_locking_mechanism_probe',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, verdict)


if __name__ == '__main__':
    main()
