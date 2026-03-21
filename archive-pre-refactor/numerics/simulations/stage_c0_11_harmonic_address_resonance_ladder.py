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
    plot_grades,
    plot_trace,
    plot_trajectories,
)
from stage23_6_dk_phase_sensitivity_braid_exchange import (
    connected_components,
    edge_grid,
    evolve,
    grade_exchange_signal,
    grade_weights,
    pair_phase_difference,
    packet_state_with_skew,
)
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions
from stage_c0_3_incidence_noise_robustness import endpoint_rewire_step
from stage_c0_4_controlled_motif_injection import (
    adjacency_connected_sets,
    choose_seed_patch,
    copy_adjacency,
    hub_specs,
    inject_balanced_motif,
    triangle_specs,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_11_harmonic_address_resonance_ladder_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_11_Harmonic_Address_Resonance_Ladder_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'run_order',
    'variant_label',
    'address_mode',
    'structural_control',
    'left_address',
    'right_address',
    'address_difference',
    'address_band_difference',
    'resolution',
    'topology_class',
    'grade_exchange_coherence',
    'resonance_alignment_index',
    'topology_stability_flag',
    'refinement_topology_class',
    'refinement_resonance_alignment_index',
    'binding_persistence',
    'flow_concentration_index',
    'transport_span',
    'channel_count',
    'loop_count',
    'phase_alignment_metric',
    'structural_edge_delta',
    'ladder_band_label',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.11 harmonic-address resonance ladder scan.')
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


def make_rng(run_id: str) -> np.random.Generator:
    seed = 5003 + sum((idx + 1) * ord(ch) for idx, ch in enumerate(run_id))
    return np.random.default_rng(seed)


def wrap_phase(angle: float) -> float:
    return float((angle + math.pi) % (2.0 * math.pi) - math.pi)


def harmonic_phase(address: int, modulus: int) -> float:
    return float((2.0 * math.pi / modulus) * address)


def raw_address_difference(a: int, b: int, modulus: int) -> int:
    return int((int(b) - int(a)) % modulus)


def minimal_address_difference(a: int, b: int, modulus: int) -> int:
    diff = abs(int(a) - int(b)) % modulus
    return int(min(diff, modulus - diff))


def operator_adjacency_lists(operator: sp.csr_matrix) -> list[list[int]]:
    adjacency = [set() for _ in range(operator.shape[0])]
    coo = operator.tocoo()
    for row, col, value in zip(coo.row, coo.col, coo.data):
        if row == col:
            continue
        if abs(value) <= 1.0e-12:
            continue
        adjacency[int(row)].add(int(col))
        adjacency[int(col)].add(int(row))
    return [sorted(items) for items in adjacency]


def median_offdiag_abs(operator: sp.csr_matrix) -> float:
    coo = operator.tocoo()
    vals = [abs(complex(value)) for row, col, value in zip(coo.row, coo.col, coo.data) if row != col and abs(value) > 1.0e-12]
    if not vals:
        return 1.0
    return float(np.median(np.asarray(vals, dtype=float)))


def combined_initial_weights(packet_states: list[np.ndarray]) -> np.ndarray:
    return np.sum(np.stack([np.abs(packet) ** 2 for packet in packet_states], axis=0), axis=0)


def adjacency_edge_set(adjacency_lists: list[list[int]]) -> set[tuple[int, int]]:
    return {
        (idx, nbr)
        for idx, nbrs in enumerate(adjacency_lists)
        for nbr in nbrs
        if nbr > idx
    }


def delta_laplacian_from_edges(
    size: int,
    baseline_edges: set[tuple[int, int]],
    modified_edges: set[tuple[int, int]],
) -> tuple[sp.csr_matrix, int]:
    added = modified_edges - baseline_edges
    removed = baseline_edges - modified_edges
    degree_delta = np.zeros(size, dtype=float)
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    for u, v in sorted(added):
        degree_delta[u] += 1.0
        degree_delta[v] += 1.0
        rows.extend([u, v])
        cols.extend([v, u])
        values.extend([-1.0, -1.0])
    for u, v in sorted(removed):
        degree_delta[u] -= 1.0
        degree_delta[v] -= 1.0
        rows.extend([u, v])
        cols.extend([v, u])
        values.extend([1.0, 1.0])
    delta = sp.coo_matrix((values, (rows, cols)), shape=(size, size), dtype=float).tocsr()
    delta = delta + sp.diags(degree_delta)
    return delta.tocsr(), len(added) + len(removed)


def apply_degree_skew(
    adjacency_lists: list[list[int]],
    avg_weights: np.ndarray,
    params: dict[str, Any],
    run_id: str,
) -> tuple[list[list[int]], int]:
    local_params = dict(params)
    local_params['hub_leaf_count'] = int(params.get('degree_skew_leaf_count', 2))
    specs = hub_specs(adjacency_lists, avg_weights, local_params)
    min_degree_floor = int(params.get('min_degree_floor', 1))
    for spec in specs:
        payload = inject_balanced_motif(adjacency_lists, spec, min_degree_floor=min_degree_floor)
        if payload is not None:
            return payload['adjacency_lists'], len(payload.get('added_edges', [])) + len(payload.get('removed_edges', []))

    rng = make_rng(run_id)
    seed, patch_nodes, _ = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=6, params=local_params)
    adj = copy_adjacency(adjacency_lists)
    operations = 0
    for _ in range(6):
        moved = endpoint_rewire_step(
            adj,
            rng,
            min_degree_floor=min_degree_floor,
            source_pool=patch_nodes,
            target_pool=[seed],
        )
        if moved is None:
            continue
        operations += 1
    modified = [sorted(items) for items in adj]
    if adjacency_connected_sets(adj) and operations > 0:
        return modified, operations
    return adjacency_lists, 0


def apply_triangle_motif(
    adjacency_lists: list[list[int]],
    avg_weights: np.ndarray,
    params: dict[str, Any],
) -> tuple[list[list[int]], int]:
    specs = triangle_specs(adjacency_lists, avg_weights, params)
    min_degree_floor = int(params.get('min_degree_floor', 1))
    for spec in specs:
        payload = inject_balanced_motif(adjacency_lists, spec, min_degree_floor=min_degree_floor)
        if payload is not None:
            return payload['adjacency_lists'], len(payload.get('added_edges', [])) + len(payload.get('removed_edges', []))
    return adjacency_lists, 0


def structural_operator(
    base_operator: sp.csr_matrix,
    packet_states: list[np.ndarray],
    structural_control: str,
    params: dict[str, Any],
    run_id: str,
) -> tuple[sp.csr_matrix, int]:
    if structural_control == 'baseline':
        return base_operator, 0

    adjacency_lists = operator_adjacency_lists(base_operator)
    avg_weights = combined_initial_weights(packet_states)
    if structural_control == 'degree_skew':
        modified_adj, edge_delta = apply_degree_skew(adjacency_lists, avg_weights, params, run_id)
    elif structural_control == 'triangle_motif':
        modified_adj, edge_delta = apply_triangle_motif(adjacency_lists, avg_weights, params)
    else:
        raise ValueError(f'unsupported structural_control: {structural_control}')

    if edge_delta <= 0:
        return base_operator, 0

    baseline_edges = adjacency_edge_set(adjacency_lists)
    modified_edges = adjacency_edge_set(modified_adj)
    delta_laplacian, edge_delta = delta_laplacian_from_edges(base_operator.shape[0], baseline_edges, modified_edges)
    scale = float(params.get('graph_delta_scale_fraction', 0.18)) * median_offdiag_abs(base_operator)
    operator = (base_operator + scale * delta_laplacian).tocsr()
    return operator, edge_delta


def resonance_alignment_index(phase_diffs: list[float], expected_phase: float) -> float:
    if not phase_diffs:
        return 0.0
    values = [0.5 * (1.0 + math.cos(wrap_phase(value - expected_phase))) for value in phase_diffs]
    return float(np.mean(np.asarray(values, dtype=float)))


def phase_alignment_metric(phase_diffs: list[float]) -> float:
    if not phase_diffs:
        return 0.0
    mean_abs = float(np.mean(np.abs(np.asarray(phase_diffs, dtype=float))))
    return float(np.clip(1.0 - mean_abs / math.pi, 0.0, 1.0))


def classify_topology(
    braid_flag: bool,
    concentration: float,
    coherence: float,
    binding_persistence: float,
    alignment: float,
) -> str:
    if braid_flag and concentration >= 0.88 and coherence >= 0.44 and alignment >= 0.68:
        return 'braid_like_exchange'
    if binding_persistence >= 0.92 and concentration >= 0.90:
        return 'localized_capture'
    if concentration >= 0.84 and coherence >= 0.36 and alignment >= 0.45:
        return 'smeared_transfer'
    return 'dispersive_pass'


def ladder_band_label(address_difference: int) -> str:
    if address_difference == 0:
        return 'band_0_exact'
    if address_difference in (1, 4):
        return 'band_1_adjacent'
    if address_difference in (2, 3):
        return 'band_2_offset'
    return 'band_random'


def resolve_addresses(run: dict[str, Any], modulus: int) -> tuple[int, int]:
    if str(run['address_mode']) != 'randomized':
        return int(run['left_address']), int(run['right_address'])
    rng = make_rng(str(run['run_id']))
    left = int(rng.integers(0, modulus))
    right = int(rng.integers(0, modulus))
    return left, right


def simulate_case(run: dict[str, Any], common: dict[str, Any], resolution: int, include_plots: bool) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    epsilon = float(common['epsilon'])
    sigma = float(common['mean_width'])
    amplitude = float(common['amplitude'])
    t_final = float(common['t_final'])
    dt_scale = float(common['dt_scale'])
    kick_cycles = float(common['kick_cycles'])
    beta = float(common['beta'])
    skew = float(common['local_phase_skew_fraction'])
    separation = float(common['separation'])
    modulus = int(common['harmonic_modulus'])

    left_address, right_address = resolve_addresses(run, modulus)
    left_phase = harmonic_phase(left_address, modulus)
    right_phase = harmonic_phase(right_address, modulus)
    expected_phase = wrap_phase(right_phase - left_phase)

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
            phase_offset=left_phase if idx == 0 else right_phase,
            kick_vector=kicks[idx],
            kick_cycles=kick_cycles,
            phase_skew_fraction=0.0 if idx == 0 else skew,
        )
        for idx in range(2)
    ]

    operator, structural_edge_delta = structural_operator(
        base_operator=complex_data.dirac_kahler,
        packet_states=packet_states,
        structural_control=str(run['structural_control']),
        params=common,
        run_id=str(run['run_id']),
    )

    psi0 = np.sum(packet_states, axis=0)
    dt = first_order_dt(complex_data.dirac_kahler, dt_scale)
    steps = max(2, int(math.ceil(t_final / dt)))
    states = evolve(operator, block_sizes, psi0, dt, steps, beta)
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
    flat = np.sort(avg_edge_grid.ravel())[::-1]
    total = float(np.sum(avg_edge_grid))
    top_k = max(1, int(math.ceil(0.1 * flat.size)))
    concentration = float(np.sum(flat[:top_k]) / max(total, 1.0e-12))

    mean_pair_distances, min_pair_distances = mean_pair_separation_series(center_histories)
    close_threshold = 2.0 * sigma
    composite_lifetime = float(np.sum(mean_pair_distances <= close_threshold) * dt)
    binding_persistence = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    final_mean_separation = float(mean_pair_distances[-1]) if mean_pair_distances.size else 0.0
    min_separation = float(np.min(min_pair_distances)) if min_pair_distances.size else 0.0
    post_collision_trend = final_mean_separation - min_separation

    braid_flag = abs(float(np.mean(braid_votes)) - 0.5) > 0.25 if braid_votes else False
    alignment = resonance_alignment_index(phase_diffs, expected_phase)
    topology_class = classify_topology(
        braid_flag=braid_flag,
        concentration=concentration,
        coherence=coherence,
        binding_persistence=binding_persistence,
        alignment=alignment,
    )
    phase_align = phase_alignment_metric(phase_diffs)
    threshold = max(0.15 * float(np.max(avg_edge_grid)), float(np.mean(avg_edge_grid) + np.std(avg_edge_grid)))
    mask = avg_edge_grid >= threshold
    components = connected_components(mask)
    loop_count = 0
    transport_span = float(max(0.0, final_mean_separation - min_separation))
    address_difference = raw_address_difference(left_address, right_address, modulus)
    band_difference = minimal_address_difference(left_address, right_address, modulus)

    row = {
        'run_id': str(run['run_id']),
        'geometry_id': 'tight_clustered_pair',
        'phase_id': str(run['variant_label']),
        'resolution': resolution,
        'graph_type': common['graph_type'],
        'packet_count': 2,
        'initial_peak_count': int(raw_peak_counts[0]) if raw_peak_counts else 0,
        'collision_label': topology_class,
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
        'new_collision_class': int(topology_class == 'braid_like_exchange'),
        'gate_met': 0,
        'promoted_followup': 0,
        'variant_label': str(run['variant_label']),
        'address_mode': str(run['address_mode']),
        'structural_control': str(run['structural_control']),
        'left_address': left_address,
        'right_address': right_address,
        'address_difference': address_difference,
        'address_band_difference': band_difference,
        'topology_class': topology_class,
        'grade_exchange_coherence': coherence,
        'resonance_alignment_index': alignment,
        'topology_stability_flag': 0,
        'refinement_topology_class': 'pending',
        'refinement_resonance_alignment_index': 0.0,
        'flow_concentration_index': concentration,
        'transport_span': transport_span,
        'channel_count': components,
        'loop_count': loop_count,
        'phase_alignment_metric': phase_align,
        'structural_edge_delta': structural_edge_delta,
        'ladder_band_label': ladder_band_label(address_difference),
        'notes': str(run['notes']),
    }

    detail = {
        'run': run,
        'metrics': row,
        'times': times,
        'phase_differences': phase_diffs,
        'grade_weights': grade_hist.tolist(),
        'grade_exchange_trace': transfer_signal.tolist(),
        'peak_tracks': [hist.tolist() for hist in center_histories],
        'flow_concentration_trace': [float(np.sum(np.sort(g.ravel())[::-1][:max(1, int(math.ceil(0.1 * g.size)))]) / max(np.sum(g), 1.0e-12)) for g in edge_grids],
        'expected_phase_offset': expected_phase,
    }

    plot_paths: list[Path] = []
    if include_plots:
        stem = f"stage_c0_11_{run['run_id']}"
        traj_path = WORK_PLOT_DIR / f'{stem}_topology_trajectory.png'
        flow_path = WORK_PLOT_DIR / f'{stem}_flow_concentration_trace.png'
        grade_path = WORK_PLOT_DIR / f'{stem}_grade_exchange_trace.png'
        grades_path = WORK_PLOT_DIR / f'{stem}_grade_weights.png'
        plot_trajectories(traj_path, str(run['run_id']), center_histories, 'tight_clustered_pair')
        plot_trace(flow_path, str(run['run_id']), times, np.asarray(detail['flow_concentration_trace'], dtype=float), 'flow concentration', 'Stage C0.11 flow concentration trace')
        plot_trace(grade_path, str(run['run_id']), times, transfer_signal, 'grade exchange', 'Stage C0.11 grade-exchange trace')
        plot_grades(grades_path, str(run['run_id']), times, grade_hist)
        plot_paths.extend([traj_path, flow_path, grade_path, grades_path])
    return row, detail, plot_paths


def plot_topology_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    mapping = {
        'braid_like_exchange': 0,
        'smeared_transfer': 1,
        'localized_capture': 2,
        'dispersive_pass': 3,
    }
    matrix = np.full((3, 6), np.nan, dtype=float)
    text: dict[tuple[int, int], str] = {}
    for row in rows:
        label = str(row['variant_label'])
        topo = str(row['topology_class'])
        if label == 'same_address_baseline':
            i, j = 0, 0
        elif label == 'difference_1':
            i, j = 0, 1
        elif label == 'difference_2':
            i, j = 0, 2
        elif label == 'difference_3':
            i, j = 0, 3
        elif label == 'difference_4':
            i, j = 0, 4
        elif label == 'random_pairing_control':
            i, j = 0, 5
        elif label == 'difference_2_degree_skew':
            i, j = 1, 2
        elif label == 'difference_3_motif':
            i, j = 1, 3
        elif label == 'same_address_refinement':
            i, j = 2, 0
        else:
            continue
        matrix[i, j] = mapping[topo]
        text[(i, j)] = topo.replace('_', ' ')
    fig, ax = plt.subplots(figsize=(10.0, 4.8))
    im = ax.imshow(matrix, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_yticks([0, 1, 2], ['baseline ladder', 'structural controls', 'refinement'])
    ax.set_xticks(range(6), ['d0', 'd1', 'd2', 'd3', 'd4', 'rand'])
    ax.set_title('Stage C0.11 topology class matrix vs address difference')
    for (i, j), label in text.items():
        ax.text(j, i, label, ha='center', va='center', color='white', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_resonance_alignment_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=lambda item: int(item['run_order']) if 'run_order' in item else 0)
    labels = [str(row['variant_label']) for row in ordered]
    alignment = [float(row['resonance_alignment_index']) for row in ordered]
    coherence = [float(row['grade_exchange_coherence']) for row in ordered]
    colors = ['tab:blue' if str(row['topology_class']) == 'braid_like_exchange' else 'tab:orange' for row in ordered]
    fig, axes = plt.subplots(2, 1, figsize=(11.0, 7.2), sharex=True)
    axes[0].bar(range(len(ordered)), alignment, color=colors)
    axes[0].set_ylabel('resonance alignment')
    axes[0].set_title('Stage C0.11 resonance alignment panel')
    axes[0].grid(axis='y', alpha=0.25)
    axes[1].bar(range(len(ordered)), coherence, color='tab:green')
    axes[1].set_ylabel('grade exchange coherence')
    axes[1].set_xticks(range(len(ordered)), labels, rotation=45, ha='right')
    axes[1].grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def ladder_summary(rows: list[dict[str, Any]]) -> tuple[str, str]:
    baseline_rows = {int(row['address_difference']): row for row in rows if str(row['structural_control']) == 'baseline' and str(row['variant_label']) in {'same_address_baseline', 'difference_1', 'difference_2', 'difference_3', 'difference_4'}}
    if len(baseline_rows) < 5:
        return 'Inconclusive', 'The baseline address-difference ladder is incomplete, so no stable resonance-band claim is justified.'

    topologies = {diff: str(row['topology_class']) for diff, row in baseline_rows.items()}
    unique_count = len(set(topologies.values()))
    cyclic_pairing = topologies[1] == topologies[4] and topologies[2] == topologies[3]
    if unique_count == 1:
        return 'No', 'The topology class stays flat across address difference, so no harmonic resonance ladder appears at the topology level.'

    transitions = [f'{left}->{right}' for left, right in zip(range(5), range(1, 5)) if topologies[left] != topologies[right]]
    if cyclic_pairing:
        alignment_14 = 0.5 * (float(baseline_rows[1]['resonance_alignment_index']) + float(baseline_rows[4]['resonance_alignment_index']))
        alignment_23 = 0.5 * (float(baseline_rows[2]['resonance_alignment_index']) + float(baseline_rows[3]['resonance_alignment_index']))
        extra = ', '.join(transitions) if transitions else 'no sharp internal transitions'
        if unique_count == 2 and topologies[1] == topologies[2]:
            return 'Yes', f'A harmonic resonance ladder is present, but it is split across observables: topology shows 2 visible bands, with d0 isolated from the dispersive d1-d4 sector, while resonance alignment further stratifies that sector into d1/d4 ({alignment_14:.3f}) versus d2/d3 ({alignment_23:.3f}). Sensitive transitions appear at {extra}.'
        band_count = len({topologies[0], topologies[1], topologies[2]})
        return 'Yes', f'A harmonic resonance ladder is present with {band_count} visible topology bands. The clean phase bands are d0, d1/d4, and d2/d3, with sensitive transitions at {extra}.'
    return 'Inconclusive', 'Topology varies across address difference, but it does not settle into a clean cyclic ladder pattern.'


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], verdict: str, note: str) -> None:
    lines = [
        '# Stage C0.11 Harmonic Address Resonance Ladder v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Ladder verdict: `{verdict}`',
        '',
        'Per-run summary:',
        '',
    ]
    ordered = sorted(rows, key=lambda item: int(item['run_order']))
    for row in ordered:
        lines.extend([
            f"- `{row['run_id']}`",
            f"  - address difference: `{row['address_difference']}`",
            f"  - structural control: `{row['structural_control']}`",
            f"  - topology class: `{row['topology_class']}`",
            f"  - resonance alignment index: `{row['resonance_alignment_index']:.4f}`",
            f"  - grade exchange coherence: `{row['grade_exchange_coherence']:.4f}`",
            f"  - topology stability flag: `{row['topology_stability_flag']}`",
            f"  - refinement topology class: `{row['refinement_topology_class']}`",
        ])
    lines.extend([
        '',
        'Classification note:',
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

    row_order = {str(run['run_id']): int(run['run_order']) for run in runsheet['runs']}

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    refinement_resolution = int(common['refinement_resolution'])

    for run in runs:
        row, detail, run_plots = simulate_case(run, common, resolution=int(common['resolution']), include_plots=True)
        row['run_order'] = row_order[str(run['run_id'])]
        refine_row, _refine_detail, _ = simulate_case(run, common, resolution=refinement_resolution, include_plots=False)
        row['topology_stability_flag'] = int(str(row['topology_class']) == str(refine_row['topology_class']))
        row['refinement_topology_class'] = str(refine_row['topology_class'])
        row['refinement_resonance_alignment_index'] = float(refine_row['resonance_alignment_index'])
        detail['metrics'] = row
        rows.append(row)
        details.append(detail)
        plot_paths.extend(run_plots)

    rows.sort(key=lambda item: int(item['run_order']))
    details.sort(key=lambda item: int(item['metrics']['run_order']))

    matrix_path = WORK_PLOT_DIR / 'stage_c0_11_topology_class_matrix_vs_address_difference.png'
    resonance_path = WORK_PLOT_DIR / 'stage_c0_11_resonance_alignment_panel.png'
    plot_topology_matrix(matrix_path, rows)
    plot_resonance_alignment_panel(resonance_path, rows)
    plot_paths.extend([matrix_path, resonance_path])

    verdict, note = ladder_summary(rows)
    summary = {
        'topology_counts': dict(Counter(str(row['topology_class']) for row in rows)),
        'stability_count': int(sum(int(row['topology_stability_flag']) for row in rows)),
        'verdict': verdict,
        'classification_note': note,
    }
    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': summary,
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage_c0_11_harmonic_address_resonance_ladder',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, verdict, note)


if __name__ == '__main__':
    main()
