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
    packet_state_with_skew,
)
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions
from stage_c0_4_controlled_motif_injection import (
    add_edge,
    adjacency_connected_sets,
    choose_seed_patch,
    copy_adjacency,
    edge_support,
    remove_edge,
    triangle_specs,
    inject_balanced_motif,
)
from stage_c0_11_harmonic_address_resonance_ladder import (
    adjacency_edge_set,
    apply_degree_skew,
    combined_initial_weights,
    delta_laplacian_from_edges,
    median_offdiag_abs,
    operator_adjacency_lists,
)
from stage_c0_12_harmonic_detuning_continuum_scan import (
    address_activity_statistics,
    address_weighted_operator,
    build_address_labels,
    classify_trace_state,
    phase_alignment_metric,
    transfer_smear_index,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_13_topology_phase_boundary_refinement_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_13_Topology_Phase_Boundary_Refinement_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'run_order',
    'family_id',
    'family_label',
    'motif_protection_strength',
    'degree_skew',
    'grade_coupling_amplitude',
    'delta',
    'right_address',
    'topology_class',
    'topology_stability_index',
    'flow_concentration_index',
    'address_selectivity_index',
    'grade_exchange_asymmetry',
    'local_channel_count',
    'recirculation_score',
    'boundary_sharpness_metric',
    'adjacent_class_variance',
    'phase_alignment_metric',
    'address_sorting_score',
    'topology_survival_time',
    'transfer_asymmetry',
    'structural_edge_delta',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.13 topology phase-boundary refinement scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def expand_runs(runsheet: dict[str, Any]) -> list[dict[str, Any]]:
    runs: list[dict[str, Any]] = []
    order = 0
    for family in sorted(runsheet['families'], key=lambda item: int(item['family_order'])):
        for delta in runsheet['detuning_samples']:
            delta_str = str(delta).replace('.', 'p')
            run_id = f"C013_{family['family_id']}_d{delta_str}"
            runs.append(
                {
                    'run_id': run_id,
                    'run_order': order,
                    'delta': float(delta),
                    **family,
                }
            )
            order += 1
    return runs


def selected_runs(runs: list[dict[str, Any]], run_ids: list[str]) -> list[dict[str, Any]]:
    if not run_ids:
        return runs
    wanted = set(run_ids)
    return [run for run in runs if str(run['run_id']) in wanted]


def topology_code(label: str) -> int:
    return {
        'address_protected_braid': 0,
        'braid_like_exchange': 1,
        'transfer_smeared': 2,
        'unresolved_mixed': 3,
    }[label]


def classify_final_topology(
    concentration: float,
    coherence: float,
    braid_flag: bool,
    selectivity: float,
    structural_edge_delta: int,
) -> str:
    if braid_flag and concentration >= 0.89 and coherence >= 0.44 and selectivity >= 0.88:
        return 'address_protected_braid'
    if structural_edge_delta > 0 and braid_flag and concentration >= 0.895 and coherence >= 0.45 and selectivity >= 0.84:
        return 'braid_like_exchange'
    if concentration >= 0.84 and coherence >= 0.34:
        return 'transfer_smeared'
    return 'unresolved_mixed'


def topology_trace_labels(
    edge_grids: list[np.ndarray],
    braid_votes: list[float],
    transfer_signal: np.ndarray,
    selectivity: float,
) -> tuple[list[str], list[float]]:
    labels: list[str] = []
    values: list[float] = []
    peak_transfer = float(np.max(transfer_signal)) if transfer_signal.size else 0.0
    for idx, grid in enumerate(edge_grids):
        total = float(np.sum(grid))
        flat = np.sort(grid.ravel())[::-1]
        top_k = max(1, int(math.ceil(0.1 * flat.size)))
        concentration = float(np.sum(flat[:top_k]) / max(total, 1.0e-12))
        braid_flag = abs(float(np.mean(braid_votes[: idx + 1])) - 0.5) > 0.25 if braid_votes[: idx + 1] else False
        transfer_norm = float(transfer_signal[idx] / peak_transfer) if peak_transfer > 1.0e-12 and idx < transfer_signal.size else 0.0
        coarse = classify_trace_state(concentration, transfer_norm, braid_flag, selectivity)
        if coarse == 'braid_like_exchange' and selectivity >= 0.88:
            label = 'address_protected_braid'
        elif coarse == 'braid_like_exchange':
            label = 'braid_like_exchange'
        elif coarse == 'transfer_smeared':
            label = 'transfer_smeared'
        else:
            label = 'unresolved_mixed'
        labels.append(label)
        values.append(float(topology_code(label)))
    return labels, values


def transfer_asymmetry(avg_edge_grid: np.ndarray) -> float:
    half = avg_edge_grid.shape[0] // 2
    left = float(np.sum(avg_edge_grid[:half, :]))
    right = float(np.sum(avg_edge_grid[half:, :]))
    total = left + right
    if total <= 1.0e-12:
        return 0.0
    return float(abs(left - right) / total)


def apply_triangle_weak_plus(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> tuple[list[list[int]], int]:
    specs = triangle_specs(adjacency_lists, avg_weights, params)
    min_degree_floor = int(params.get('min_degree_floor', 1))
    for spec in specs:
        payload = inject_balanced_motif(adjacency_lists, spec, min_degree_floor=min_degree_floor)
        if payload is not None:
            edge_delta = len(payload.get('added_edges', [])) + len(payload.get('removed_edges', []))
            return payload['adjacency_lists'], edge_delta
    return adjacency_lists, 0


def apply_triangle_weak_minus(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> tuple[list[list[int]], int]:
    min_degree_floor = int(params.get('min_degree_floor', 1))
    seed, patch_nodes, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=4, params=params)
    patch_set = set(int(node) for node in patch_nodes)
    adj_sets = copy_adjacency(adjacency_lists)

    candidate_edges: list[tuple[float, int, int]] = []
    for u in patch_nodes:
        for v in adj_sets[u]:
            if v <= u or v not in patch_set:
                continue
            if len(adj_sets[u]) <= min_degree_floor or len(adj_sets[v]) <= min_degree_floor:
                continue
            support = edge_support(adj_sets, u, v)
            if support <= 0:
                continue
            score = float(support) + float(avg_weights[u]) + float(avg_weights[v])
            candidate_edges.append((score, int(u), int(v)))
    candidate_edges.sort(reverse=True)

    outside_nodes = [idx for idx in range(len(adjacency_lists)) if idx not in patch_set]
    for _score, u, v in candidate_edges:
        remove_edge(adj_sets, u, v)
        if not adjacency_connected_sets(adj_sets):
            add_edge(adj_sets, u, v)
            continue
        added = None
        for a_idx, a in enumerate(outside_nodes):
            for b in outside_nodes[a_idx + 1:]:
                if b in adj_sets[a]:
                    continue
                add_edge(adj_sets, a, b)
                added = (a, b)
                break
            if added is not None:
                break
        if added is None:
            add_edge(adj_sets, u, v)
            continue
        return [sorted(items) for items in adj_sets], 2
    return adjacency_lists, 0


def structural_operator(
    base_operator: sp.csr_matrix,
    packet_states: list[np.ndarray],
    motif_strength: str,
    degree_skew: str,
    params: dict[str, Any],
    run_id: str,
) -> tuple[sp.csr_matrix, int]:
    adjacency_lists = operator_adjacency_lists(base_operator)
    avg_weights = combined_initial_weights(packet_states)
    modified_adj = adjacency_lists
    total_edge_delta = 0

    if motif_strength == 'weak_plus':
        modified_adj, edge_delta = apply_triangle_weak_plus(modified_adj, avg_weights, params)
        total_edge_delta += edge_delta
    elif motif_strength == 'weak_minus':
        modified_adj, edge_delta = apply_triangle_weak_minus(modified_adj, avg_weights, params)
        total_edge_delta += edge_delta

    if degree_skew == 'mild':
        local = dict(params)
        local['degree_skew_leaf_count'] = int(params.get('degree_skew_leaf_count', 2))
        modified_adj, edge_delta = apply_degree_skew(modified_adj, avg_weights, local, run_id)
        total_edge_delta += edge_delta
    elif degree_skew == 'asymmetric':
        local = dict(params)
        local['degree_skew_leaf_count'] = int(params.get('degree_skew_leaf_count_asymmetric', 4))
        modified_adj, edge_delta = apply_degree_skew(modified_adj, avg_weights, local, run_id + '_asym')
        total_edge_delta += edge_delta

    if total_edge_delta <= 0:
        return base_operator, 0

    baseline_edges = adjacency_edge_set(adjacency_lists)
    modified_edges = adjacency_edge_set(modified_adj)
    delta_laplacian, _ = delta_laplacian_from_edges(base_operator.shape[0], baseline_edges, modified_edges)
    scale = float(params.get('graph_delta_scale_fraction', 0.18)) * median_offdiag_abs(base_operator)
    operator = (base_operator + scale * delta_laplacian).tocsr()
    return operator, total_edge_delta


def local_channel_count(avg_edge_grid: np.ndarray) -> int:
    threshold = max(0.15 * float(np.max(avg_edge_grid)), float(np.mean(avg_edge_grid) + np.std(avg_edge_grid)))
    mask = avg_edge_grid >= threshold
    return connected_components(mask)


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
    beta = float(run['beta'])
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
    structured_operator, structural_edge_delta = structural_operator(
        base_operator=base_operator,
        packet_states=packet_states,
        motif_strength=str(run['motif_protection_strength']),
        degree_skew=str(run['degree_skew']),
        params=common,
        run_id=str(run['run_id']),
    )
    effective_operator = address_weighted_operator(structured_operator, labels, modulus, eta_match, eta_mismatch)
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
    _same_activity, sorting_score, selectivity = address_activity_statistics(
        states=states,
        operator=effective_operator,
        labels=labels,
        modulus=modulus,
        eta_match=eta_match,
        eta_mismatch=eta_mismatch,
    )
    topology_class = classify_final_topology(concentration, coherence, braid_flag, selectivity, structural_edge_delta)
    trace_labels, trace_values = topology_trace_labels(edge_grids, braid_votes, transfer_signal, selectivity)
    stability_index = float(np.mean([label == topology_class for label in trace_labels])) if trace_labels else 0.0
    topology_survival_time = float(sum(1 for label in trace_labels if label in {'address_protected_braid', 'braid_like_exchange'}) * dt)
    asymmetry = transfer_asymmetry(avg_edge_grid)
    phase_align = phase_alignment_metric(phase_diffs)
    smear_index = transfer_smear_index(concentration, coherence)
    grade_asymmetry = float(abs(grade_hist[-1, 0] - grade_hist[-1, 1]))

    row = {
        'run_id': str(run['run_id']),
        'geometry_id': 'tight_clustered_pair',
        'phase_id': f"delta_{delta:.3f}",
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
        'new_collision_class': int(topology_class in {'address_protected_braid', 'braid_like_exchange'}),
        'gate_met': 0,
        'promoted_followup': 0,
        'run_order': int(run['run_order']),
        'family_id': str(run['family_id']),
        'family_label': str(run['family_label']),
        'motif_protection_strength': str(run['motif_protection_strength']),
        'degree_skew': str(run['degree_skew']),
        'grade_coupling_amplitude': str(run['grade_coupling_amplitude']),
        'delta': delta,
        'right_address': right_address,
        'topology_class': topology_class,
        'topology_stability_index': stability_index,
        'flow_concentration_index': concentration,
        'address_selectivity_index': selectivity,
        'grade_exchange_asymmetry': grade_asymmetry,
        'local_channel_count': local_channel_count(avg_edge_grid),
        'recirculation_score': 0.0,
        'boundary_sharpness_metric': 0.0,
        'adjacent_class_variance': 0.0,
        'phase_alignment_metric': phase_align,
        'address_sorting_score': sorting_score,
        'topology_survival_time': topology_survival_time,
        'transfer_asymmetry': asymmetry,
        'structural_edge_delta': structural_edge_delta,
        'notes': str(run['notes']),
    }

    detail = {
        'run': run,
        'metrics': row,
        'times': times,
        'topology_trace': trace_labels,
        'topology_trace_values': trace_values,
        'grade_exchange_trace': transfer_signal.tolist(),
        'flow_concentration_trace': [float(np.sum(np.sort(g.ravel())[::-1][:max(1, int(math.ceil(0.1 * g.size)))]) / max(np.sum(g), 1.0e-12)) for g in edge_grids],
    }

    plot_paths: list[Path] = []
    topology_trace_path = WORK_PLOT_DIR / f"stage_c0_13_{run['run_id']}_topology_vs_time_trace.png"
    selectivity_trace_path = WORK_PLOT_DIR / f"stage_c0_13_{run['run_id']}_selectivity_trace.png"
    plot_trace(topology_trace_path, str(run['run_id']), times[: len(trace_values)], np.asarray(trace_values, dtype=float), 'topology code', 'Stage C0.13 topology-vs-time trace')
    selectivity_trace = np.linspace(1.0, selectivity, num=len(times))
    plot_trace(selectivity_trace_path, str(run['run_id']), times, selectivity_trace, 'address selectivity', 'Stage C0.13 selectivity-vs-time trace')
    plot_paths.extend([topology_trace_path, selectivity_trace_path])
    return row, detail, plot_paths


def compute_family_diagnostics(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row['family_id'])].append(row)

    diagnostics: dict[str, dict[str, Any]] = {}
    for family_id, family_rows in grouped.items():
        family_rows.sort(key=lambda item: float(item['delta']))
        codes = np.asarray([topology_code(str(row['topology_class'])) for row in family_rows], dtype=float)
        deltas = np.asarray([float(row['delta']) for row in family_rows], dtype=float)
        selectivity = np.asarray([float(row['address_selectivity_index']) for row in family_rows], dtype=float)
        stability = np.asarray([float(row['topology_stability_index']) for row in family_rows], dtype=float)

        local_sharpness = np.zeros(len(family_rows), dtype=float)
        local_variance = np.zeros(len(family_rows), dtype=float)
        for idx in range(len(family_rows)):
            neigh_codes = [codes[idx]]
            if idx > 0:
                delta_step = max(deltas[idx] - deltas[idx - 1], 1.0e-9)
                sharp = abs(codes[idx] - codes[idx - 1]) + abs(selectivity[idx] - selectivity[idx - 1]) + abs(stability[idx] - stability[idx - 1])
                local_sharpness[idx] += sharp / delta_step
                neigh_codes.append(codes[idx - 1])
            if idx + 1 < len(family_rows):
                delta_step = max(deltas[idx + 1] - deltas[idx], 1.0e-9)
                sharp = abs(codes[idx + 1] - codes[idx]) + abs(selectivity[idx + 1] - selectivity[idx]) + abs(stability[idx + 1] - stability[idx])
                local_sharpness[idx] += sharp / delta_step
                neigh_codes.append(codes[idx + 1])
            local_variance[idx] = float(np.var(np.asarray(neigh_codes, dtype=float)))
            family_rows[idx]['boundary_sharpness_metric'] = float(local_sharpness[idx])
            family_rows[idx]['adjacent_class_variance'] = float(local_variance[idx])

        resolved = False
        boundary_delta: float | None = None
        for start in range(len(family_rows) - 2):
            band = codes[start:start + 3]
            if not np.all(band == band[0]):
                continue
            left_ok = start == 0 or codes[start - 1] != band[0]
            right_ok = start + 3 >= len(codes) or codes[start + 3] != band[0]
            if left_ok and right_ok:
                resolved = True
                if start + 3 < len(codes):
                    boundary_delta = float(deltas[start + 3])
                elif start > 0:
                    boundary_delta = float(deltas[start])
                break

        diagnostics[family_id] = {
            'resolved': resolved,
            'boundary_delta': boundary_delta,
            'max_sharpness': float(np.max(local_sharpness)) if len(local_sharpness) else 0.0,
            'mean_adjacent_variance': float(np.mean(local_variance)) if len(local_variance) else 0.0,
            'rows': family_rows,
        }
    return diagnostics


def plot_topology_phase_diagram(path: Path, rows: list[dict[str, Any]]) -> None:
    family_order = sorted({(str(row['family_id']), str(row['family_label'])) for row in rows}, key=lambda item: next(int(r['run_order']) for r in rows if str(r['family_id']) == item[0]))
    family_ids = [item[0] for item in family_order]
    family_labels = [item[1] for item in family_order]
    deltas = sorted({float(row['delta']) for row in rows})
    matrix = np.full((len(family_ids), len(deltas)), np.nan, dtype=float)
    text: dict[tuple[int, int], str] = {}
    for i, family_id in enumerate(family_ids):
        for j, delta in enumerate(deltas):
            row = next(item for item in rows if str(item['family_id']) == family_id and abs(float(item['delta']) - delta) < 1.0e-9)
            matrix[i, j] = float(topology_code(str(row['topology_class'])))
            text[(i, j)] = str(row['topology_class']).replace('_', ' ')
    fig, ax = plt.subplots(figsize=(11.0, 5.8))
    im = ax.imshow(matrix, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_yticks(range(len(family_labels)), family_labels)
    ax.set_xticks(range(len(deltas)), [f'{delta:.2f}' for delta in deltas])
    ax.set_xlabel('harmonic detuning')
    ax.set_title('Stage C0.13 topology phase diagram')
    for (i, j), label in text.items():
        ax.text(j, i, label, ha='center', va='center', color='white', fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_boundary_sharpness_heatmap(path: Path, diagnostics: dict[str, dict[str, Any]]) -> None:
    family_ids = list(diagnostics.keys())
    family_labels = [str(diagnostics[family_id]['rows'][0]['family_label']) for family_id in family_ids]
    deltas = [float(row['delta']) for row in diagnostics[family_ids[0]]['rows']]
    matrix = np.asarray([[float(row['boundary_sharpness_metric']) for row in diagnostics[family_id]['rows']] for family_id in family_ids], dtype=float)
    fig, ax = plt.subplots(figsize=(11.0, 5.8))
    im = ax.imshow(matrix, aspect='auto', cmap='magma')
    ax.set_yticks(range(len(family_labels)), family_labels)
    ax.set_xticks(range(len(deltas)), [f'{delta:.2f}' for delta in deltas])
    ax.set_xlabel('harmonic detuning')
    ax.set_title('Stage C0.13 boundary-sharpness heatmap')
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_topology_trace_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row['family_id'])].append(row)
    fig, ax = plt.subplots(figsize=(11.0, 5.4))
    mapping = {
        'address_protected_braid': 0,
        'braid_like_exchange': 1,
        'transfer_smeared': 2,
        'unresolved_mixed': 3,
    }
    for family_id, family_rows in grouped.items():
        family_rows.sort(key=lambda item: float(item['delta']))
        xs = [float(row['delta']) for row in family_rows]
        ys = [mapping[str(row['topology_class'])] for row in family_rows]
        label = str(family_rows[0]['family_label'])
        ax.plot(xs, ys, marker='o', label=label)
    ax.set_yticks([0, 1, 2, 3], ['prot braid', 'braid', 'smear', 'mixed'])
    ax.set_xlabel('harmonic detuning')
    ax.set_title('Stage C0.13 topology vs detuning trace')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc='best')
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_selectivity_trace_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row['family_id'])].append(row)
    fig, ax = plt.subplots(figsize=(11.0, 5.4))
    for family_id, family_rows in grouped.items():
        family_rows.sort(key=lambda item: float(item['delta']))
        xs = [float(row['delta']) for row in family_rows]
        ys = [float(row['address_selectivity_index']) for row in family_rows]
        label = str(family_rows[0]['family_label'])
        ax.plot(xs, ys, marker='o', label=label)
    ax.set_xlabel('harmonic detuning')
    ax.set_ylabel('address selectivity index')
    ax.set_title('Stage C0.13 selectivity vs detuning trace')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc='best')
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], diagnostics: dict[str, dict[str, Any]], verdict: str) -> None:
    lines = [
        '# Stage C0.13 Topology Phase Boundary Refinement v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Global read: `{verdict}`',
        '',
        '## Protection-sensitivity table',
    ]
    for family_id, info in diagnostics.items():
        row0 = info['rows'][0]
        lines.append(
            f"- `{row0['family_label']}`: resolved=`{int(info['resolved'])}`, boundary_delta=`{info['boundary_delta']}`, max_sharpness=`{info['max_sharpness']:.3f}`, adjacent_variance=`{info['mean_adjacent_variance']:.3f}`"
        )
    lines.extend([
        '',
        '## Per-run summary',
    ])
    for row in sorted(rows, key=lambda item: int(item['run_order'])):
        lines.extend([
            f"- `{row['run_id']}`",
            f"  - family: `{row['family_label']}`",
            f"  - delta: `{row['delta']:.3f}`",
            f"  - topology class: `{row['topology_class']}`",
            f"  - topology stability index: `{row['topology_stability_index']:.4f}`",
            f"  - address selectivity index: `{row['address_selectivity_index']:.4f}`",
            f"  - boundary sharpness metric: `{row['boundary_sharpness_metric']:.4f}`",
        ])
    lines.extend([
        '',
        '## Plots',
    ])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')
    NOTE_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def global_verdict(diagnostics: dict[str, dict[str, Any]]) -> str:
    resolved_count = sum(int(info['resolved']) for info in diagnostics.values())
    boundary_count = len({info['boundary_delta'] for info in diagnostics.values() if info['boundary_delta'] is not None})
    if resolved_count >= 4 and boundary_count <= 2:
        return 'harmonic addressing produces a sharp topology phase transition'
    if resolved_count >= 2:
        return 'harmonic addressing produces a finite crossover corridor'
    return 'harmonic addressing produces fragmented multi-window protection'


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    runs = selected_runs(expand_runs(runsheet), args.run_ids)
    common = runsheet['common_fields']

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        row, detail, run_plots = simulate_run(run, common)
        rows.append(row)
        details.append(detail)
        plot_paths.extend(run_plots)

    rows.sort(key=lambda item: int(item['run_order']))
    diagnostics = compute_family_diagnostics(rows)
    for detail in details:
        run_id = str(detail['run']['run_id'])
        detail['metrics'] = next(row for row in rows if str(row['run_id']) == run_id)
    details.sort(key=lambda item: int(item['metrics']['run_order']))

    phase_path = WORK_PLOT_DIR / 'stage_c0_13_topology_phase_diagram_panel.png'
    sharpness_path = WORK_PLOT_DIR / 'stage_c0_13_boundary_sharpness_heatmap.png'
    topology_trace_path = WORK_PLOT_DIR / 'stage_c0_13_topology_vs_detuning_trace.png'
    selectivity_trace_path = WORK_PLOT_DIR / 'stage_c0_13_selectivity_vs_detuning_trace.png'
    plot_topology_phase_diagram(phase_path, rows)
    plot_boundary_sharpness_heatmap(sharpness_path, diagnostics)
    plot_topology_trace_panel(topology_trace_path, rows)
    plot_selectivity_trace_panel(selectivity_trace_path, rows)
    plot_paths.extend([phase_path, sharpness_path, topology_trace_path, selectivity_trace_path])

    verdict = global_verdict(diagnostics)
    summary = {
        'topology_counts': dict(Counter(str(row['topology_class']) for row in rows)),
        'family_diagnostics': {
            family_id: {
                'resolved': bool(info['resolved']),
                'boundary_delta': info['boundary_delta'],
                'max_sharpness': info['max_sharpness'],
                'mean_adjacent_variance': info['mean_adjacent_variance'],
            }
            for family_id, info in diagnostics.items()
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
        'stage_c0_13_topology_phase_boundary_refinement',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, diagnostics, verdict)


if __name__ == '__main__':
    main()
