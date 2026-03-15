#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import itertools
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, load_stage10_defaults, plt, save_atlas_payload
from stage_c0_1_combinatorial_kernel_probe import build_line_graph_cache
from stage_c0_3_incidence_noise_robustness import build_incidence_perturbation, compute_shells_from_adjacency
from stage_c0_4_controlled_motif_injection import choose_seed_patch, edge_key, inject_balanced_motif
from stage_c0_5_emergent_compression_preorder import (
    adjacency_sets,
    analyze_support,
    apply_pass_separation,
    canonical_partition,
    focused_partition,
    initial_edge_weights,
    local_square_counts,
    local_triangle_counts,
    multi_source_shell_distance,
    pairwise_partition_mismatch,
    partition_map,
    quantize_metric,
    restricted_neighbors,
    select_active_support,
    top_seed_edges,
    survivor_partition,
    survivor_support,
)
from stage_c0_5_emergent_sequencing import count_seed_components, load_reference_rows, topology_label_for_support
from stage_c0_6_incidence_seeded_mismatch_accumulation import forced_mismatch_ledger, survivor_overlay_plot

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_8_weaker_kernel_bi_root_motif_injection_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_8_weaker_kernel_bi_root_motif_injection'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

NOISE_IDS = [
    'ultra_weak_rewiring',
    'weak_rewiring',
    'moderate_rewiring',
    'moderate_deletion_insertion',
    'clustered_rewiring_patch',
    'random_localized_defect_pair',
    'degree_biased_perturbation',
    'connectivity_edge_stress_test',
]
MOTIF_IDS = [
    'bi_root_triangle_cluster',
    'bi_root_diamond',
    'bi_root_hub_micro_star',
    'two_layer_stacked_motif',
]

CSV_FIELDS = [
    'run_id',
    'representative_id',
    'forcing_protocol_id',
    'forcing_protocol_label',
    'protected_topology_label',
    'reference_interaction_label',
    'reference_coarse_label',
    'resolution',
    'weaker_kernel_id',
    'weaker_kernel_hop_cap',
    'initial_active_support',
    'initial_seed_count',
    'steps_completed',
    'tau_sequencing_depth',
    'mismatch_activation_step',
    'protected_arrow_step',
    'protected_arrow_step_distribution',
    'protected_arrow_step_under_weaker_inconsistency',
    'persistent_step',
    'final_topology_label',
    'final_active_persistence_class_count',
    'final_seed_component_count',
    'final_survivor_count',
    'obstruction_rank_final',
    'kernel_update_consistency_failure_count',
    'kernel_update_consistency_failure_max',
    'kernel_update_consistency_failure_initial',
    'kernel_update_consistency_failure_final',
    'positive_mismatch_density',
    'max_mismatch_density',
    'final_mismatch_density',
    'positive_mismatch_steps',
    'bi_root_motif_correlation',
    'protected_blocked_merge_total',
    'static_collapse_flag',
    'sequencing_class',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.8 weaker-kernel inconsistency forcing scan.')
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


def lookup_by_id(items: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    return {str(item[key]): item for item in items}


def forcing_protocol_order(protocol_id: str) -> int:
    return {
        'noise_seeded_weaker_inconsistency': 0,
        'bi_root_motif_seeded_weaker_inconsistency': 1,
        'hybrid_seeded_weaker_inconsistency': 2,
    }.get(protocol_id, 99)


def payload_touched_nodes(payload: dict[str, Any], forcing_kind: str) -> list[int]:
    if forcing_kind == 'noise':
        touch_counts = np.asarray(payload.get('touch_counts', []), dtype=int)
        return [int(idx) for idx in np.flatnonzero(touch_counts > 0)]
    touched: set[int] = set(int(value) for value in payload.get('motif_nodes', []))
    for key in ('added_edges', 'removed_edges', 'desired_edges'):
        for edge in payload.get(key, []):
            if len(edge) == 2:
                touched.add(int(edge[0]))
                touched.add(int(edge[1]))
    return sorted(touched)


def specs_missing_edges(adjacency_lists: list[list[int]], desired: set[tuple[int, int]]) -> list[tuple[int, int]]:
    return sorted(edge for edge in desired if edge[1] not in adjacency_lists[edge[0]])


def bi_root_triangle_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=5, params=params)
    subset = candidates[: min(len(candidates), 8)]
    specs: list[tuple[float, dict[str, Any]]] = []
    for roots in itertools.combinations(subset[: min(len(subset), 6)], 2):
        remaining = [int(node) for node in subset if int(node) not in roots]
        for branches in itertools.combinations(remaining, 2):
            a, b = [int(value) for value in roots]
            c, d = [int(value) for value in branches]
            desired = {
                edge_key(a, b),
                edge_key(a, c),
                edge_key(b, c),
                edge_key(a, d),
                edge_key(b, d),
            }
            missing = specs_missing_edges(adjacency_lists, desired)
            if not missing or len(missing) > 4:
                continue
            motif_nodes = sorted({a, b, c, d})
            score = sum(float(avg_weights[idx]) for idx in motif_nodes) - 0.12 * len(missing)
            specs.append(
                (
                    score,
                    {
                        'motif_nodes': motif_nodes,
                        'desired_edges': sorted(desired),
                        'seed_nodes': [a, b],
                        'patch_nodes': patch,
                        'missing_edges': missing,
                    },
                )
            )
    specs.sort(key=lambda item: item[0], reverse=True)
    return [item[1] for item in specs[:24]]


def bi_root_diamond_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=4, params=params)
    subset = candidates[: min(len(candidates), 8)]
    specs: list[tuple[float, dict[str, Any]]] = []
    for roots in itertools.combinations(subset[: min(len(subset), 6)], 2):
        remaining = [int(node) for node in subset if int(node) not in roots]
        for branches in itertools.combinations(remaining, 2):
            a, b = [int(value) for value in roots]
            c, d = [int(value) for value in branches]
            desired = {
                edge_key(a, c),
                edge_key(b, c),
                edge_key(a, d),
                edge_key(b, d),
                edge_key(c, d),
            }
            missing = specs_missing_edges(adjacency_lists, desired)
            if not missing or len(missing) > 4:
                continue
            motif_nodes = sorted({a, b, c, d})
            score = sum(float(avg_weights[idx]) for idx in motif_nodes) - 0.10 * len(missing)
            specs.append(
                (
                    score,
                    {
                        'motif_nodes': motif_nodes,
                        'desired_edges': sorted(desired),
                        'seed_nodes': [a, b],
                        'patch_nodes': patch,
                        'missing_edges': missing,
                    },
                )
            )
    specs.sort(key=lambda item: item[0], reverse=True)
    return [item[1] for item in specs[:24]]


def bi_root_hub_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=6, params=params)
    subset = [int(node) for node in candidates[: min(len(candidates), 8)]]
    specs: list[tuple[float, dict[str, Any]]] = []
    for roots in itertools.combinations(subset[: min(len(subset), 6)], 2):
        remaining = [node for node in subset if node not in roots]
        if len(remaining) < 4:
            continue
        a, b = [int(value) for value in roots]
        leaves = remaining[:4]
        desired = {
            edge_key(a, b),
            edge_key(a, leaves[0]),
            edge_key(a, leaves[1]),
            edge_key(b, leaves[2]),
            edge_key(b, leaves[3]),
        }
        missing = specs_missing_edges(adjacency_lists, desired)
        if not missing or len(missing) > 4:
            continue
        motif_nodes = sorted({a, b, *leaves})
        score = sum(float(avg_weights[idx]) for idx in motif_nodes) - 0.10 * len(missing)
        specs.append(
            (
                score,
                {
                    'motif_nodes': motif_nodes,
                    'desired_edges': sorted(desired),
                    'seed_nodes': [a, b],
                    'patch_nodes': patch,
                    'missing_edges': missing,
                },
            )
        )
    specs.sort(key=lambda item: item[0], reverse=True)
    return [item[1] for item in specs[:24]]


def stacked_two_layer_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=5, params=params)
    subset = [int(node) for node in candidates[: min(len(candidates), 8)]]
    specs: list[tuple[float, dict[str, Any]]] = []
    for combo in itertools.combinations(subset, 5):
        a, b, c, d, e = [int(value) for value in combo]
        desired = {
            edge_key(a, b),
            edge_key(a, c),
            edge_key(b, c),
            edge_key(c, d),
            edge_key(c, e),
            edge_key(d, e),
        }
        missing = specs_missing_edges(adjacency_lists, desired)
        if not missing or len(missing) > 5:
            continue
        motif_nodes = sorted({a, b, c, d, e})
        score = sum(float(avg_weights[idx]) for idx in motif_nodes) - 0.08 * len(missing)
        specs.append(
            (
                score,
                {
                    'motif_nodes': motif_nodes,
                    'desired_edges': sorted(desired),
                    'seed_nodes': [a, b, c],
                    'patch_nodes': patch,
                    'missing_edges': missing,
                },
            )
        )
    specs.sort(key=lambda item: item[0], reverse=True)
    return [item[1] for item in specs[:24]]


def build_bi_root_motif_injection(
    adjacency_lists: list[list[int]],
    motif_id: str,
    avg_weights: np.ndarray,
    params: dict[str, Any],
) -> dict[str, Any]:
    if motif_id == 'bi_root_triangle_cluster':
        specs = bi_root_triangle_specs(adjacency_lists, avg_weights, params)
    elif motif_id == 'bi_root_diamond':
        specs = bi_root_diamond_specs(adjacency_lists, avg_weights, params)
    elif motif_id == 'bi_root_hub_micro_star':
        specs = bi_root_hub_specs(adjacency_lists, avg_weights, params)
    elif motif_id == 'two_layer_stacked_motif':
        specs = stacked_two_layer_specs(adjacency_lists, avg_weights, params)
    else:
        raise ValueError(f'unsupported motif id: {motif_id}')

    min_degree_floor = int(params.get('min_degree_floor', 3))
    for spec in specs:
        payload = inject_balanced_motif(adjacency_lists, spec, min_degree_floor=min_degree_floor)
        if payload is not None:
            return payload
    raise RuntimeError(f'failed to inject motif {motif_id}')


def motif_payload_from_adjacency(
    adjacency_lists: list[list[int]],
    motif_id: str,
    avg_weights: np.ndarray,
    params: dict[str, Any],
    forcing_kind: str,
    forcing_item_id: str,
    extra_touched_nodes: list[int] | None = None,
) -> dict[str, Any]:
    payload = dict(build_bi_root_motif_injection(adjacency_lists, motif_id, avg_weights, params))
    touched = set(payload_touched_nodes(payload, 'motif'))
    touched.update(int(node) for node in (extra_touched_nodes or []))
    payload['forcing_kind'] = forcing_kind
    payload['forcing_item_id'] = str(forcing_item_id)
    payload['touched_nodes'] = sorted(touched)
    payload['shells'] = compute_shells_from_adjacency(payload['adjacency_lists'], max_shell=int(params.get('max_shell', 4)))
    return payload


def build_noise_payloads(data: Any, params: dict[str, Any]) -> dict[str, dict[str, Any]]:
    payloads: dict[str, dict[str, Any]] = {}
    for perturbation_id in NOISE_IDS:
        payload = dict(build_incidence_perturbation(data, perturbation_id, params))
        payload['forcing_kind'] = 'noise'
        payload['forcing_item_id'] = str(perturbation_id)
        payload['touched_nodes'] = payload_touched_nodes(payload, 'noise')
        payload['shells'] = compute_shells_from_adjacency(payload['adjacency_lists'], max_shell=int(params.get('max_shell', 4)))
        payloads[str(perturbation_id)] = payload
    return payloads


def build_baseline_motif_payloads(
    baseline_adj_lists: list[list[int]],
    avg_weights: np.ndarray,
    params: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    payloads: dict[str, dict[str, Any]] = {}
    for motif_id in MOTIF_IDS:
        payloads[str(motif_id)] = motif_payload_from_adjacency(
            baseline_adj_lists,
            motif_id,
            avg_weights,
            params,
            forcing_kind='motif',
            forcing_item_id=str(motif_id),
        )
    return payloads


def forcing_candidates_for_step(
    forcing_protocol_id: str,
    step_idx: int,
    noise_payloads: dict[str, dict[str, Any]],
    motif_payloads: dict[str, dict[str, Any]],
    avg_weights: np.ndarray,
    params: dict[str, Any],
) -> list[dict[str, Any]]:
    noise_id = NOISE_IDS[min(step_idx, len(NOISE_IDS) - 1)]
    motif_id = MOTIF_IDS[step_idx % len(MOTIF_IDS)]
    if forcing_protocol_id == 'noise_seeded_weaker_inconsistency':
        return [noise_payloads[noise_id]]
    if forcing_protocol_id == 'bi_root_motif_seeded_weaker_inconsistency':
        return [motif_payloads[motif_id]]
    if forcing_protocol_id == 'hybrid_seeded_weaker_inconsistency':
        noise_payload = noise_payloads[noise_id]
        combined = motif_payload_from_adjacency(
            noise_payload['adjacency_lists'],
            motif_id,
            avg_weights,
            params,
            forcing_kind='hybrid',
            forcing_item_id=f'{noise_id}+{motif_id}',
            extra_touched_nodes=list(noise_payload['touched_nodes']),
        )
        combined['source_noise_id'] = str(noise_id)
        combined['source_motif_id'] = str(motif_id)
        return [combined]
    raise ValueError(f'unsupported forcing protocol id: {forcing_protocol_id}')


def analyze_support_c08(
    support: list[int],
    weights: np.ndarray,
    full_adj_sets: list[set[int]],
    full_shells: list[dict[int, int]],
    seed_shell_distance: np.ndarray,
    params: dict[str, Any],
) -> dict[str, Any]:
    base_params = dict(params)
    base_params['max_shell'] = min(int(params.get('max_shell', 4)), 3)
    state = dict(analyze_support(support, weights, full_adj_sets, full_shells, seed_shell_distance, base_params))
    support = [int(node) for node in state.get('support', [])]
    empty_partition: tuple[tuple[int, ...], ...] = tuple()
    if not support:
        state.update(
            {
                'bounded_path_partition': empty_partition,
                'bounded_path_map': {},
                'bounded_path_survivor_count': 0,
                'bounded_path_equivalence_count': 0,
            }
        )
        return state

    max_hop = int(params.get('bounded_path_max_hop', 4))
    levels = int(params.get('bounded_path_levels', 8))
    restricted = restricted_neighbors(support, full_adj_sets)
    triangle_counts = local_triangle_counts(restricted)
    square_counts = local_square_counts(restricted)
    degrees = {int(node): len(restricted[int(node)]) for node in support}
    support_set = set(support)

    path_mass: dict[int, float] = {}
    path_degree_mass: dict[int, float] = {}
    path_reach: dict[int, int] = {}
    shell_counts: dict[int, list[int]] = {}

    for node in support:
        shell_profile = [0 for _ in range(max_hop)]
        weighted_mass = 0.0
        weighted_degree = 0.0
        for nbr, dist in full_shells[int(node)].items():
            if nbr not in support_set or dist < 1 or dist > max_hop:
                continue
            factor = 1.0 / float(dist + 1)
            weighted_mass += factor * float(weights[int(nbr)])
            weighted_degree += factor * float(degrees.get(int(nbr), 0))
            shell_profile[dist - 1] += 1
        path_mass[int(node)] = weighted_mass
        path_degree_mass[int(node)] = weighted_degree
        path_reach[int(node)] = int(sum(shell_profile))
        shell_counts[int(node)] = shell_profile

    max_weight = float(np.max(weights[np.asarray(support, dtype=int)]))
    max_degree = float(max(degrees.values())) if degrees else 0.0
    max_triangle = float(max(triangle_counts.values())) if triangle_counts else 0.0
    max_square = float(max(square_counts.values())) if square_counts else 0.0
    max_path_mass = max(path_mass.values()) if path_mass else 0.0
    max_path_degree = max(path_degree_mass.values()) if path_degree_mass else 0.0
    max_reach = float(max(path_reach.values())) if path_reach else 0.0
    shell_max = [float(max(shell_counts[node][idx] for node in support)) for idx in range(max_hop)]

    groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for node in support:
        seed_bin = int(min(int(seed_shell_distance[int(node)]), max_hop + 1))
        shell_bins = tuple(
            quantize_metric(float(shell_counts[int(node)][idx]), shell_max[idx], levels=4) for idx in range(max_hop)
        )
        signature = (
            seed_bin,
            quantize_metric(float(weights[int(node)]), max_weight, levels=levels),
            quantize_metric(float(degrees[int(node)]), max_degree, levels=levels),
            quantize_metric(float(triangle_counts[int(node)]), max_triangle, levels=4),
            quantize_metric(float(square_counts[int(node)]), max_square, levels=4),
            quantize_metric(float(path_mass[int(node)]), max_path_mass, levels=levels),
            quantize_metric(float(path_degree_mass[int(node)]), max_path_degree, levels=levels),
            quantize_metric(float(path_reach[int(node)]), max_reach, levels=levels),
        ) + shell_bins
        groups[signature].append(int(node))

    bounded_partition = canonical_partition(list(groups.values()))
    bounded_map = partition_map(bounded_partition)
    state.update(
        {
            'bounded_path_partition': bounded_partition,
            'bounded_path_map': bounded_map,
            'bounded_path_survivor_count': len(survivor_support(bounded_partition)),
            'bounded_path_equivalence_count': sum(1 for cls in bounded_partition if len(cls) > 1),
        }
    )
    return state


def next_partition_for_kernel_c08(state: dict[str, Any], kernel_name: str) -> tuple[tuple[int, ...], ...]:
    if kernel_name == 'adjacency':
        return state['adjacency_partition']
    if kernel_name == 'bounded_path':
        return state['bounded_path_partition']
    if kernel_name == 'graph_shell':
        return state['graph_shell_partition']
    raise ValueError(f'unsupported kernel name: {kernel_name}')


def class_map_for_kernel_c08(state: dict[str, Any], kernel_name: str) -> dict[int, tuple[int, ...]]:
    if kernel_name == 'adjacency':
        return state['adjacency_map']
    if kernel_name == 'bounded_path':
        return state['bounded_path_map']
    if kernel_name == 'graph_shell':
        return state['graph_shell_map']
    raise ValueError(f'unsupported kernel name: {kernel_name}')


def evaluate_candidate_c08(
    state: dict[str, Any],
    current_partition: tuple[tuple[int, ...], ...],
    focus_nodes: list[int],
    kernel_name: str,
    pass_active: bool,
    weights: np.ndarray,
    full_adj_sets: list[set[int]],
    full_shells: list[dict[int, int]],
    seed_shell_distance: np.ndarray,
    params: dict[str, Any],
) -> dict[str, Any]:
    raw_partition = next_partition_for_kernel_c08(state, kernel_name)
    blocked_merges = 0
    final_partition = raw_partition
    if pass_active:
        final_partition, blocked_merges = apply_pass_separation(raw_partition, current_partition)
    final_partition = focused_partition(final_partition, focus_nodes)
    next_support = survivor_support(final_partition)
    next_state = analyze_support_c08(next_support, weights, full_adj_sets, full_shells, seed_shell_distance, params)
    return {
        'kernel_name': kernel_name,
        'raw_partition': raw_partition,
        'final_partition': survivor_partition(final_partition),
        'next_support': next_support,
        'next_state': next_state,
        'blocked_merges': int(blocked_merges),
    }


def restricted_partition(partition: tuple[tuple[int, ...], ...], nodes: set[int]) -> tuple[tuple[int, ...], ...]:
    return canonical_partition(
        [
            tuple(sorted(int(node) for node in cls if int(node) in nodes))
            for cls in partition
            if any(int(node) in nodes for node in cls)
        ]
    )


def partition_inconsistency_count(
    left_partition: tuple[tuple[int, ...], ...],
    right_partition: tuple[tuple[int, ...], ...],
) -> int:
    left_nodes = set(int(node) for cls in left_partition for node in cls)
    right_nodes = set(int(node) for cls in right_partition for node in cls)
    common_nodes = sorted(left_nodes.intersection(right_nodes))
    support_delta = len(left_nodes.symmetric_difference(right_nodes))
    if len(common_nodes) < 2:
        return int(support_delta)
    left_map = partition_map(restricted_partition(left_partition, set(common_nodes)))
    right_map = partition_map(restricted_partition(right_partition, set(common_nodes)))
    pairwise_count, _contrib = pairwise_partition_mismatch(common_nodes, left_map, right_map)
    return int(pairwise_count + support_delta)


def follow_up_candidate(
    first_candidate: dict[str, Any],
    next_kernel: str,
    seed_edges: list[int],
    weights: np.ndarray,
    full_adj_sets: list[set[int]],
    full_shells: list[dict[int, int]],
    seed_shell_distance: np.ndarray,
    params: dict[str, Any],
) -> dict[str, Any]:
    focus_nodes = [int(node) for node in seed_edges if int(node) in first_candidate['next_support']] or list(first_candidate['next_support'])
    return evaluate_candidate_c08(
        first_candidate['next_state'],
        first_candidate['final_partition'],
        focus_nodes,
        next_kernel,
        True,
        weights,
        full_adj_sets,
        full_shells,
        seed_shell_distance,
        params,
    )


def motif_ledger_correlation(node_count: int, motif_nodes: list[int], ledger: list[int]) -> float:
    if not motif_nodes or not ledger:
        return 0.0
    motif_indicator = np.zeros(node_count, dtype=float)
    ledger_indicator = np.zeros(node_count, dtype=float)
    motif_indicator[np.asarray(sorted(set(int(node) for node in motif_nodes)), dtype=int)] = 1.0
    ledger_indicator[np.asarray(sorted(set(int(node) for node in ledger)), dtype=int)] = 1.0
    if float(np.std(motif_indicator)) <= 1.0e-18 or float(np.std(ledger_indicator)) <= 1.0e-18:
        return 0.0
    return float(np.corrcoef(motif_indicator, ledger_indicator)[0, 1])


def evaluate_payload(
    payload: dict[str, Any],
    current_support: list[int],
    current_partition: tuple[tuple[int, ...], ...],
    current_state: dict[str, Any],
    current_label: str,
    current_kernel: str,
    seed_edges: list[int],
    weights: np.ndarray,
    params: dict[str, Any],
    protected_label: str,
) -> list[dict[str, Any]]:
    adjacency_lists = payload['adjacency_lists']
    full_adj_sets = adjacency_sets(adjacency_lists)
    full_shells = payload['shells']
    seed_shell_distance = multi_source_shell_distance(adjacency_lists, seed_edges, int(params.get('seed_shell_cap', 4)))
    forced_state = analyze_support_c08(current_support, weights, full_adj_sets, full_shells, seed_shell_distance, params)
    focus_nodes = [int(node) for node in seed_edges if int(node) in current_support] or list(current_support)

    adjacency_candidate = evaluate_candidate_c08(
        forced_state,
        current_partition,
        focus_nodes,
        'adjacency',
        True,
        weights,
        full_adj_sets,
        full_shells,
        seed_shell_distance,
        params,
    )
    bounded_candidate = evaluate_candidate_c08(
        forced_state,
        current_partition,
        focus_nodes,
        'bounded_path',
        True,
        weights,
        full_adj_sets,
        full_shells,
        seed_shell_distance,
        params,
    )
    adjacency_label = topology_label_for_support(
        adjacency_candidate['next_support'],
        adjacency_candidate['final_partition'],
        adjacency_candidate['next_state'],
        seed_edges,
        full_adj_sets,
        protected_label,
    )
    bounded_label = topology_label_for_support(
        bounded_candidate['next_support'],
        bounded_candidate['final_partition'],
        bounded_candidate['next_state'],
        seed_edges,
        full_adj_sets,
        protected_label,
    )

    adjacency_then_bounded = follow_up_candidate(
        adjacency_candidate,
        'bounded_path',
        seed_edges,
        weights,
        full_adj_sets,
        full_shells,
        seed_shell_distance,
        params,
    )
    bounded_then_adjacency = follow_up_candidate(
        bounded_candidate,
        'adjacency',
        seed_edges,
        weights,
        full_adj_sets,
        full_shells,
        seed_shell_distance,
        params,
    )
    failure_count = partition_inconsistency_count(
        adjacency_then_bounded['final_partition'],
        bounded_then_adjacency['final_partition'],
    )

    results: list[dict[str, Any]] = []
    for kernel_name, candidate, label in (
        ('adjacency', adjacency_candidate, adjacency_label),
        ('bounded_path', bounded_candidate, bounded_label),
    ):
        ledger = forced_mismatch_ledger(
            prev_state=current_state,
            prev_support=current_support,
            next_state=candidate['next_state'],
            next_support=candidate['next_support'],
            touched_nodes=list(payload['touched_nodes']),
        )
        density = float(len(ledger) / max(len(adjacency_lists), 1))
        motif_corr = motif_ledger_correlation(len(adjacency_lists), list(payload.get('motif_nodes', [])), ledger)
        results.append(
            {
                'forcing_kind': str(payload['forcing_kind']),
                'forcing_item_id': str(payload['forcing_item_id']),
                'current_topology_label': str(current_label),
                'adjacency_topology_label': str(adjacency_label),
                'bounded_path_topology_label': str(bounded_label),
                'active_persistence_class_count': int(len({current_label, adjacency_label, bounded_label})),
                'kernel_name': str(kernel_name),
                'candidate': candidate,
                'chosen_topology_label': str(label),
                'preserves_protected_label': int(label == protected_label),
                'mismatch_ledger': ledger,
                'mismatch_density': density,
                'bi_root_motif_correlation': float(motif_corr),
                'kernel_update_consistency_failure_count': int(failure_count),
                'forced_state': forced_state,
                'touched_nodes': list(payload['touched_nodes']),
                'motif_nodes': [int(node) for node in payload.get('motif_nodes', [])],
                'graph_info': {
                    'adjacency_lists': adjacency_lists,
                    'full_adj_sets': full_adj_sets,
                    'full_shells': full_shells,
                    'seed_shell_distance': seed_shell_distance,
                },
                'current_kernel_match': int(kernel_name == current_kernel),
            }
        )
    return results


def blocked_hold_candidate(
    template: dict[str, Any],
    current_support: list[int],
    current_partition: tuple[tuple[int, ...], ...],
    current_state: dict[str, Any],
    current_label: str,
    current_kernel: str,
    current_graph_info: dict[str, Any],
) -> dict[str, Any]:
    ledger = forced_mismatch_ledger(
        prev_state=current_state,
        prev_support=current_support,
        next_state=template['forced_state'],
        next_support=current_support,
        touched_nodes=template['touched_nodes'],
    )
    density = float(len(ledger) / max(len(template['graph_info']['adjacency_lists']), 1))
    motif_corr = motif_ledger_correlation(
        len(template['graph_info']['adjacency_lists']),
        list(template.get('motif_nodes', [])),
        ledger,
    )
    return {
        'forcing_kind': str(template['forcing_kind']),
        'forcing_item_id': str(template['forcing_item_id']),
        'current_topology_label': str(current_label),
        'adjacency_topology_label': str(template['adjacency_topology_label']),
        'bounded_path_topology_label': str(template['bounded_path_topology_label']),
        'active_persistence_class_count': int(template['active_persistence_class_count']),
        'kernel_name': str(current_kernel),
        'candidate': {
            'blocked_merges': 1,
            'next_support': list(current_support),
            'next_state': current_state,
            'final_partition': current_partition,
        },
        'chosen_topology_label': str(current_label),
        'preserves_protected_label': 1,
        'mismatch_ledger': ledger,
        'mismatch_density': density,
        'bi_root_motif_correlation': float(motif_corr),
        'kernel_update_consistency_failure_count': int(template['kernel_update_consistency_failure_count']),
        'forced_state': template['forced_state'],
        'touched_nodes': list(template['touched_nodes']),
        'motif_nodes': list(template.get('motif_nodes', [])),
        'graph_info': current_graph_info,
        'current_kernel_match': 1,
    }


def select_candidate(
    candidates: list[dict[str, Any]],
    current_support: list[int],
    current_partition: tuple[tuple[int, ...], ...],
    current_state: dict[str, Any],
    current_label: str,
    current_kernel: str,
    current_graph_info: dict[str, Any],
    inconsistency_active: bool,
) -> dict[str, Any]:
    pool = candidates
    if inconsistency_active:
        target_kernel = 'adjacency' if current_kernel == 'bounded_path' else 'bounded_path'
        pool = [item for item in pool if str(item['kernel_name']) == target_kernel]
    preserving = [item for item in pool if int(item['preserves_protected_label']) == 1]
    if not preserving:
        collapsed: list[dict[str, Any]] = []
        seen: set[tuple[str, str]] = set()
        for item in pool:
            key = (str(item['forcing_kind']), str(item['forcing_item_id']))
            if key in seen:
                continue
            seen.add(key)
            collapsed.append(
                blocked_hold_candidate(
                    template=item,
                    current_support=current_support,
                    current_partition=current_partition,
                    current_state=current_state,
                    current_label=current_label,
                    current_kernel=current_kernel,
                    current_graph_info=current_graph_info,
                )
            )
        pool = collapsed
    else:
        pool = preserving

    if inconsistency_active:
        ranked = sorted(
            pool,
            key=lambda item: (
                -int(item['kernel_update_consistency_failure_count']),
                -int(item['candidate']['blocked_merges'] > 0),
                -float(item['mismatch_density']),
                -abs(float(item['bi_root_motif_correlation'])),
                int(item['candidate']['next_state']['support_size']),
                -int(item['current_kernel_match']),
            ),
        )
    else:
        ranked = sorted(
            pool,
            key=lambda item: (
                -int(item['candidate']['blocked_merges'] > 0),
                -float(item['mismatch_density']),
                -abs(float(item['bi_root_motif_correlation'])),
                int(item['candidate']['next_state']['support_size']),
                -int(item['current_kernel_match']),
            ),
        )
    return ranked[0]


def sequencing_depth(trace: list[dict[str, Any]]) -> int:
    counts = [int(item['active_persistence_class_count']) for item in trace]
    for idx in range(len(counts)):
        if len(set(counts[idx:])) == 1:
            return int(trace[idx]['step_index'])
    return int(trace[-1]['step_index']) if trace else 0


def sequencing_class(summary: dict[str, Any]) -> str:
    if int(summary['static_collapse_flag']) == 1 and int(summary['protected_arrow_step_under_weaker_inconsistency']) < 0:
        return 'static_regime'
    if (
        int(summary['kernel_update_consistency_failure_count']) > 0
        and int(summary['protected_arrow_step_under_weaker_inconsistency']) >= 0
        and int(summary['persistent_step']) >= 0
    ):
        return 'weaker_inconsistency_ordered_persistent'
    if int(summary['kernel_update_consistency_failure_count']) > 0 and int(summary['protected_arrow_step_under_weaker_inconsistency']) >= 0:
        return 'weaker_inconsistency_ordered_transient'
    if int(summary['kernel_update_consistency_failure_count']) > 0:
        return 'weaker_inconsistency_stable_without_arrow'
    if float(summary['positive_mismatch_density']) > 0.0 and int(summary['protected_arrow_step']) >= 0:
        return 'ordered_without_weaker_inconsistency'
    if float(summary['positive_mismatch_density']) > 0.0:
        return 'rigid_positive_mismatch'
    return 'stable_without_ordering'


def weaker_trace_plot(path: Path, run_id: str, trace: list[dict[str, Any]]) -> None:
    steps = [int(item['step_index']) for item in trace]
    mismatch_density = [float(item['mismatch_density']) for item in trace]
    failure_count = [int(item['kernel_update_consistency_failure_count']) for item in trace]
    blocked = [int(item['protected_blocked_merges']) for item in trace]
    motif_corr = [float(item['bi_root_motif_correlation']) for item in trace]

    fig, axes = plt.subplots(1, 4, figsize=(17.2, 4.2))
    axes[0].plot(steps, [int(item['active_persistence_class_count']) for item in trace], marker='o', color='tab:blue')
    axes[0].set_title('Active persistence classes')
    axes[0].set_xlabel('forced step')
    axes[0].grid(alpha=0.25)

    axes[1].plot(steps, mismatch_density, marker='o', color='tab:red', label='mismatch density')
    axes[1].plot(steps, failure_count, marker='s', color='tab:purple', label='consistency failures')
    axes[1].set_title('Mismatch and failures')
    axes[1].set_xlabel('forced step')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)

    axes[2].bar(steps, blocked, alpha=0.45, color='tab:gray', label='protected blocks')
    axes[2].plot(steps, [1 if int(item['inconsistency_active']) == 1 else 0 for item in trace], marker='o', color='tab:green', label='inconsistency active')
    axes[2].set_title('Protected arrow forcing')
    axes[2].set_xlabel('forced step')
    axes[2].grid(alpha=0.25)
    axes[2].legend(fontsize=8)

    axes[3].plot(steps, motif_corr, marker='o', color='tab:orange')
    axes[3].set_title('Bi-root motif correlation')
    axes[3].set_xlabel('forced step')
    axes[3].grid(alpha=0.25)

    fig.suptitle(run_id)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def simulate_run(
    run: dict[str, Any],
    representative: dict[str, Any],
    protocol: dict[str, Any],
    runsheet: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
    protected_labels: dict[str, str],
    c01_ref: dict[tuple[str, str], dict[str, str]],
) -> tuple[dict[str, Any], list[Path]]:
    params = dict(runsheet.get('parameter_policy', {}))
    resolution = int(run['resolution'])
    boundary_type = str(run['boundary_type'])
    protected_label = str(protected_labels.get(str(run['representative_id']), representative.get('protected_topology_label', 'smeared_transfer')))

    data, weights = initial_edge_weights(representative, resolution, boundary_type, defaults, base)
    baseline_graph = build_line_graph_cache(data, max_shell=int(params.get('max_shell', 4)))
    baseline_adj_lists = baseline_graph['adjacency_lists']
    baseline_adj_sets = adjacency_sets(baseline_adj_lists)
    baseline_shells = baseline_graph['shells']

    noise_payloads = build_noise_payloads(data, params)
    avg_weights = np.asarray(weights, dtype=float)
    avg_weights = avg_weights / max(float(np.sum(avg_weights)), 1.0e-12)
    motif_payloads = build_baseline_motif_payloads(baseline_adj_lists, avg_weights, params)

    initial_support = select_active_support(weights, float(params.get('active_threshold_sigma', 1.0)), int(params.get('min_active_count', 48)))
    seed_edges = top_seed_edges(weights, initial_support, int(params.get('seed_edge_count', 12)))
    seed_shell_distance = multi_source_shell_distance(baseline_adj_lists, seed_edges, int(params.get('seed_shell_cap', 4)))
    initial_state = analyze_support_c08(initial_support, weights, baseline_adj_sets, baseline_shells, seed_shell_distance, params)

    current_partition = survivor_partition(next_partition_for_kernel_c08(initial_state, 'bounded_path'))
    current_support = survivor_support(current_partition)
    current_state = analyze_support_c08(current_support, weights, baseline_adj_sets, baseline_shells, seed_shell_distance, params)
    current_kernel = 'bounded_path'
    current_label = topology_label_for_support(current_support, current_partition, current_state, seed_edges, baseline_adj_sets, protected_label)
    current_adj_sets = baseline_adj_sets
    current_graph_info = {
        'adjacency_lists': baseline_adj_lists,
        'full_adj_sets': baseline_adj_sets,
        'full_shells': baseline_shells,
        'seed_shell_distance': seed_shell_distance,
    }

    mismatch_threshold = float(params.get('mismatch_density_threshold', 0.001))
    rho_window = int(params.get('rho_stability_window', 2))
    max_steps = int(params.get('max_steps', 8))

    trace: list[dict[str, Any]] = []
    arrow_steps: list[int] = []
    arrow_under_inconsistency: int | None = None
    mismatch_activation_step: int | None = None
    positive_mismatch_steps = 0
    positive_mismatch_density = 0.0
    max_mismatch_density = 0.0
    failure_total = 0
    failure_max = 0
    persistent_step: int | None = None
    last_cycle_rank: int | None = None
    stable_cycle_count = 0
    blocked_total = 0
    inconsistency_active = False

    for step_idx in range(max_steps):
        options = forcing_candidates_for_step(
            str(run['forcing_protocol_id']),
            step_idx,
            noise_payloads,
            motif_payloads,
            avg_weights,
            params,
        )
        candidate_records: list[dict[str, Any]] = []
        for payload in options:
            candidate_records.extend(
                evaluate_payload(
                    payload=payload,
                    current_support=current_support,
                    current_partition=current_partition,
                    current_state=current_state,
                    current_label=current_label,
                    current_kernel=current_kernel,
                    seed_edges=seed_edges,
                    weights=weights,
                    params=params,
                    protected_label=protected_label,
                )
            )
        chosen = select_candidate(
            candidate_records,
            current_support=current_support,
            current_partition=current_partition,
            current_state=current_state,
            current_label=current_label,
            current_kernel=current_kernel,
            current_graph_info=current_graph_info,
            inconsistency_active=inconsistency_active,
        )

        mismatch_density = float(chosen['mismatch_density'])
        max_mismatch_density = max(max_mismatch_density, mismatch_density)
        if mismatch_density > mismatch_threshold:
            positive_mismatch_steps += 1
            positive_mismatch_density = max(positive_mismatch_density, mismatch_density)
            if mismatch_activation_step is None:
                mismatch_activation_step = int(step_idx)

        failure_count = int(chosen['kernel_update_consistency_failure_count']) if inconsistency_active else 0
        failure_total += failure_count
        failure_max = max(failure_max, failure_count)

        blocked_merges = int(chosen['candidate']['blocked_merges'])
        if blocked_merges > 0:
            arrow_steps.append(int(step_idx + 1))
            if inconsistency_active and failure_count > 0 and arrow_under_inconsistency is None:
                arrow_under_inconsistency = int(step_idx + 1)
        blocked_total += blocked_merges

        cycle_rank = int(chosen['candidate']['next_state']['cycle_rank'])
        if last_cycle_rank == cycle_rank:
            stable_cycle_count += 1
        else:
            stable_cycle_count = 1
            last_cycle_rank = cycle_rank
        if persistent_step is None and stable_cycle_count >= rho_window:
            persistent_step = int(step_idx + 1)

        trace.append(
            {
                'step_index': int(step_idx),
                'forcing_kind': str(chosen['forcing_kind']),
                'forcing_item_id': str(chosen['forcing_item_id']),
                'current_topology_label': str(chosen['current_topology_label']),
                'adjacency_topology_label': str(chosen['adjacency_topology_label']),
                'bounded_path_topology_label': str(chosen['bounded_path_topology_label']),
                'active_persistence_class_count': int(chosen['active_persistence_class_count']),
                'chosen_kernel': str(chosen['kernel_name']),
                'chosen_topology_label': str(chosen['chosen_topology_label']),
                'mismatch_density': mismatch_density,
                'bi_root_motif_correlation': float(chosen['bi_root_motif_correlation']),
                'kernel_update_consistency_failure_count': failure_count,
                'protected_blocked_merges': blocked_merges,
                'obstruction_rank_after': cycle_rank,
                'survivor_count_after': int(len(chosen['candidate']['next_support'])),
                'seed_component_count_after': int(count_seed_components(chosen['candidate']['next_support'], seed_edges, chosen['graph_info']['full_adj_sets'])),
                'inconsistency_active': int(inconsistency_active),
            }
        )

        current_support = list(chosen['candidate']['next_support'])
        current_partition = tuple(tuple(int(node) for node in cls) for cls in chosen['candidate']['final_partition'])
        current_state = chosen['candidate']['next_state']
        current_label = str(chosen['chosen_topology_label'])
        current_kernel = str(chosen['kernel_name'])
        current_adj_sets = chosen['graph_info']['full_adj_sets']
        current_graph_info = chosen['graph_info']

        if not inconsistency_active and mismatch_activation_step is not None and step_idx >= mismatch_activation_step:
            inconsistency_active = True

    tau_depth = sequencing_depth(trace)
    reference = c01_ref.get((str(run['representative_id']), 'graph_shell'), {})
    final_seed_components = int(count_seed_components(current_support, seed_edges, current_adj_sets))
    summary = {
        'run_id': str(run['run_id']),
        'representative_id': str(run['representative_id']),
        'forcing_protocol_id': str(run['forcing_protocol_id']),
        'forcing_protocol_label': str(protocol['label']),
        'protected_topology_label': protected_label,
        'reference_interaction_label': str(reference.get('interaction_label', 'unknown')),
        'reference_coarse_label': str(reference.get('coarse_label', 'unknown')),
        'resolution': resolution,
        'weaker_kernel_id': 'bounded_path',
        'weaker_kernel_hop_cap': int(params.get('bounded_path_max_hop', 4)),
        'initial_active_support': int(len(initial_support)),
        'initial_seed_count': int(len(seed_edges)),
        'steps_completed': int(len(trace)),
        'tau_sequencing_depth': int(tau_depth),
        'mismatch_activation_step': int(mismatch_activation_step) if mismatch_activation_step is not None else -1,
        'protected_arrow_step': int(arrow_steps[0]) if arrow_steps else -1,
        'protected_arrow_step_distribution': json.dumps(arrow_steps),
        'protected_arrow_step_under_weaker_inconsistency': int(arrow_under_inconsistency) if arrow_under_inconsistency is not None else -1,
        'persistent_step': int(persistent_step) if persistent_step is not None else -1,
        'final_topology_label': str(current_label),
        'final_active_persistence_class_count': int(trace[-1]['active_persistence_class_count']) if trace else 0,
        'final_seed_component_count': final_seed_components,
        'final_survivor_count': int(len(current_support)),
        'obstruction_rank_final': int(current_state['cycle_rank']),
        'kernel_update_consistency_failure_count': int(failure_total),
        'kernel_update_consistency_failure_max': int(failure_max),
        'kernel_update_consistency_failure_initial': int(trace[0]['kernel_update_consistency_failure_count']) if trace else 0,
        'kernel_update_consistency_failure_final': int(trace[-1]['kernel_update_consistency_failure_count']) if trace else 0,
        'positive_mismatch_density': float(positive_mismatch_density),
        'max_mismatch_density': float(max_mismatch_density),
        'final_mismatch_density': float(trace[-1]['mismatch_density']) if trace else 0.0,
        'positive_mismatch_steps': int(positive_mismatch_steps),
        'bi_root_motif_correlation': float(max((abs(float(item['bi_root_motif_correlation'])) for item in trace), default=0.0)),
        'protected_blocked_merge_total': int(blocked_total),
        'static_collapse_flag': int(positive_mismatch_density <= mismatch_threshold and (arrow_under_inconsistency is None)),
        'sequencing_class': '',
        'notes': str(run['notes']),
    }
    summary['sequencing_class'] = sequencing_class(summary)

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_weaker_trace.png"
    weaker_trace_plot(trace_path, str(run['run_id']), trace)
    plot_paths.append(trace_path)

    overlay_path = WORK_PLOT_DIR / f"{run['run_id']}_support_overlay.png"
    survivor_overlay_plot(
        overlay_path,
        str(run['run_id']),
        np.asarray(data.midpoints, dtype=float),
        np.asarray(weights, dtype=float),
        initial_support,
        current_support,
        seed_edges,
        boundary_type,
    )
    plot_paths.append(overlay_path)

    result = {
        'run': run,
        'forcing_protocol': protocol,
        'summary': summary,
        'trace': trace,
        'seed_edges': [int(node) for node in seed_edges],
        'initial_active_support': [int(node) for node in initial_support],
        'final_support': [int(node) for node in current_support],
    }
    return result, plot_paths


def create_summary_plots(rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    reps = sorted({str(row['representative_id']) for row in rows})
    protocols = sorted({str(row['forcing_protocol_id']) for row in rows}, key=forcing_protocol_order)
    row_map = {(str(row['representative_id']), str(row['forcing_protocol_id'])): row for row in rows}

    class_scores = {
        'static_regime': 0,
        'stable_without_ordering': 1,
        'rigid_positive_mismatch': 2,
        'ordered_without_weaker_inconsistency': 3,
        'weaker_inconsistency_stable_without_arrow': 4,
        'weaker_inconsistency_ordered_transient': 5,
        'weaker_inconsistency_ordered_persistent': 6,
    }

    matrix_path = WORK_PLOT_DIR / 'stage_c0_8_weaker_class_matrix.png'
    matrix = np.zeros((len(reps), len(protocols)), dtype=float)
    for i, rep in enumerate(reps):
        for j, protocol in enumerate(protocols):
            matrix[i, j] = float(class_scores.get(str(row_map[(rep, protocol)]['sequencing_class']), 0))
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    im = ax.imshow(matrix, cmap='viridis', vmin=0, vmax=max(class_scores.values()))
    ax.set_xticks(range(len(protocols)))
    ax.set_xticklabels(protocols, rotation=20, ha='right')
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('C0.8 weaker-kernel class matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(matrix_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(matrix_path)

    inconsistency_path = WORK_PLOT_DIR / 'stage_c0_8_weaker_inconsistency_panel.png'
    x = np.arange(len(reps))
    width = 0.22
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6))
    for idx, protocol in enumerate(protocols):
        fail_values = [int(row_map[(rep, protocol)]['kernel_update_consistency_failure_count']) for rep in reps]
        fail_max = [int(row_map[(rep, protocol)]['kernel_update_consistency_failure_max']) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, fail_values, width=width, label=protocol)
        axes[1].bar(x + (idx - 1) * width, fail_max, width=width, label=protocol)
    axes[0].set_title('Consistency failure total')
    axes[1].set_title('Consistency failure max')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
    axes[0].legend(fontsize=8)
    fig.savefig(inconsistency_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(inconsistency_path)

    arrow_path = WORK_PLOT_DIR / 'stage_c0_8_arrow_and_motif_panel.png'
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 4.6))
    for idx, protocol in enumerate(protocols):
        arrow_values = [max(int(row_map[(rep, protocol)]['protected_arrow_step_under_weaker_inconsistency']), 0) for rep in reps]
        mismatch_values = [float(row_map[(rep, protocol)]['positive_mismatch_density']) for rep in reps]
        motif_values = [float(row_map[(rep, protocol)]['bi_root_motif_correlation']) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, arrow_values, width=width, label=protocol)
        axes[1].bar(x + (idx - 1) * width, mismatch_values, width=width, label=protocol)
        axes[2].bar(x + (idx - 1) * width, motif_values, width=width, label=protocol)
    axes[0].set_title('Arrow step under weaker inconsistency')
    axes[1].set_title('Positive mismatch density')
    axes[2].set_title('Bi-root motif correlation')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
    axes[0].legend(fontsize=8)
    fig.savefig(arrow_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(arrow_path)
    return plot_paths


def write_note(path: Path, json_rel: str, csv_rel: str, rows: list[dict[str, Any]], runsheet: dict[str, Any]) -> None:
    positive_failure_runs = sum(1 for row in rows if int(row['kernel_update_consistency_failure_count']) > 0)
    arrow_runs = sum(1 for row in rows if int(row['protected_arrow_step_under_weaker_inconsistency']) >= 0)
    lines = [
        '# Stage C0.8 Weaker Kernel Bi-Root Motif Injection v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Architecture notice: {runsheet["architecture_notice"]}',
        f'Baseline notice: {runsheet["baseline_notice"]}',
        '',
        f'Stage summary: explicit consistency failures in `{positive_failure_runs}/{len(rows)}` runs; protected arrow under weaker inconsistency in `{arrow_runs}/{len(rows)}` runs.',
        '',
        'Per-run summary:',
        '',
    ]
    for row in rows:
        lines.extend(
            [
                f"- `{row['forcing_protocol_label']}` / `{row['representative_id']}`",
                f"  - protected topology label: `{row['protected_topology_label']}`",
                f"  - sequencing class: `{row['sequencing_class']}`",
                f"  - final topology label: `{row['final_topology_label']}`",
                f"  - consistency failure total: `{row['kernel_update_consistency_failure_count']}`",
                f"  - consistency failure max: `{row['kernel_update_consistency_failure_max']}`",
                f"  - protected arrow step under weaker inconsistency: `{row['protected_arrow_step_under_weaker_inconsistency']}`",
                f"  - positive mismatch density: `{float(row['positive_mismatch_density']):.4f}`",
                f"  - bi-root motif correlation: `{float(row['bi_root_motif_correlation']):.4f}`",
            ]
        )
    lines.extend(
        [
            '',
            'Interpretation boundary:',
            '- this scan asks whether a softer bounded-path combinatorial kernel can open explicit inconsistency under the same protected sequencing rule that stayed rigid in C0.7',
            '- topology labels remain persistence-class observables and the protected topology-stability rule stays active',
            '- if consistency failures become positive and the protected arrow still fires broadly, the branch supports an irreversibility mechanism tied directly to adjacency versus weaker-kernel inconsistency',
        ]
    )
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    defaults = load_stage10_defaults()
    base = dict(runsheet['base_seed_reference'])
    representatives = lookup_by_id(runsheet['representatives'], 'representative_id')
    protocols = lookup_by_id(runsheet['forcing_protocols'], 'forcing_protocol_id')
    protected_labels, c01_ref = load_reference_rows()

    results: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in selected_runs(runsheet['runs'], args.run_ids):
        representative = representatives[str(run['representative_id'])]
        protocol = protocols[str(run['forcing_protocol_id'])]
        result, run_plots = simulate_run(run, representative, protocol, runsheet, defaults, base, protected_labels, c01_ref)
        results.append(result)
        csv_rows.append(dict(result['summary']))
        plot_paths.extend(run_plots)

    csv_rows.sort(key=lambda row: (str(row['representative_id']), forcing_protocol_order(str(row['forcing_protocol_id']))))
    plot_paths.extend(create_summary_plots(csv_rows))

    payload = {
        'stage': str(runsheet['stage']),
        'description': str(runsheet['description']),
        'architecture_notice': str(runsheet['architecture_notice']),
        'baseline_notice': str(runsheet['baseline_notice']),
        'results': results,
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug=str(runsheet.get('experiment_slug', 'stage_c0_8_weaker_kernel_bi_root_motif_injection')),
        result=payload,
        csv_rows=csv_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    note_path = ATLAS_NOTES / str(runsheet['note_name'])
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), csv_rows, runsheet)
    payload['plots'] = stamped_plots
    json_path.write_text(json.dumps(payload, indent=2), encoding='utf-8')
    print(json.dumps({'json': str(json_path), 'csv': str(csv_path), 'note': str(note_path)}, indent=2))


if __name__ == '__main__':
    main()
