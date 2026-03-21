#!/usr/bin/env python3

from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, load_stage10_defaults, plt, save_atlas_payload
from stage_c0_4_controlled_motif_injection import choose_seed_patch, edge_key, inject_balanced_motif
from stage_c0_8_weaker_kernel_bi_root_motif_injection import (
    NOISE_IDS,
    adjacency_sets,
    analyze_support_c08,
    build_line_graph_cache,
    build_noise_payloads,
    count_seed_components,
    evaluate_candidate_c08,
    forced_mismatch_ledger,
    initial_edge_weights,
    load_reference_rows,
    multi_source_shell_distance,
    partition_inconsistency_count,
    payload_touched_nodes,
    select_active_support,
    survivor_overlay_plot,
    survivor_partition,
    survivor_support,
    top_seed_edges,
    topology_label_for_support,
)
from stage_c0_3_incidence_noise_robustness import compute_shells_from_adjacency

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_9_tri_root_motif_injection_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_9_tri_root_motif_injection'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

MOTIF_IDS = [
    'tri_root_triangle_cluster',
    'tri_root_diamond',
    'tri_root_hub_micro_star',
    'three_layer_stacked_motif',
    'hybrid_bi_root_tri_root_cross_motif',
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
    'protected_arrow_step_under_tri_root_inconsistency',
    'persistent_step',
    'final_topology_label',
    'final_active_persistence_class_count',
    'final_seed_component_count',
    'final_survivor_count',
    'obstruction_rank_final',
    'kernel_update_consistency_failure_count',
    'kernel_update_consistency_failure_count_under_tri_root',
    'kernel_update_consistency_failure_max',
    'kernel_update_consistency_failure_initial',
    'kernel_update_consistency_failure_final',
    'positive_mismatch_density',
    'max_mismatch_density',
    'final_mismatch_density',
    'positive_mismatch_steps',
    'tri_root_motif_correlation',
    'irreversibility_scaling_factor',
    'protected_blocked_merge_total',
    'static_collapse_flag',
    'sequencing_class',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.9 tri-root motif escalation scan.')
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
        'noise_seeded_tri_root_inconsistency': 0,
        'tri_root_motif_seeded_inconsistency': 1,
        'hybrid_seeded_tri_root_inconsistency': 2,
    }.get(protocol_id, 99)


def specs_missing_edges(adjacency_lists: list[list[int]], desired: set[tuple[int, int]]) -> list[tuple[int, int]]:
    return sorted(edge for edge in desired if edge[1] not in adjacency_lists[edge[0]])


def tri_root_triangle_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    _seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=5, params=params)
    subset = [int(node) for node in candidates[: min(len(candidates), 8)]]
    specs: list[tuple[float, dict[str, Any]]] = []
    for roots in itertools.combinations(subset[: min(len(subset), 6)], 3):
        remaining = [node for node in subset if node not in roots]
        for branch in remaining[: min(len(remaining), 4)]:
            a, b, c = [int(value) for value in roots]
            d = int(branch)
            desired = {
                edge_key(a, b),
                edge_key(a, c),
                edge_key(b, c),
                edge_key(a, d),
                edge_key(b, d),
                edge_key(c, d),
            }
            missing = specs_missing_edges(adjacency_lists, desired)
            if not missing or len(missing) > 5:
                continue
            motif_nodes = sorted({a, b, c, d})
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


def tri_root_diamond_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    _seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=6, params=params)
    subset = [int(node) for node in candidates[: min(len(candidates), 8)]]
    specs: list[tuple[float, dict[str, Any]]] = []
    for roots in itertools.combinations(subset[: min(len(subset), 6)], 3):
        remaining = [node for node in subset if node not in roots]
        for branches in itertools.combinations(remaining[: min(len(remaining), 5)], 2):
            a, b, c = [int(value) for value in roots]
            d, e = [int(value) for value in branches]
            desired = {
                edge_key(a, d),
                edge_key(b, d),
                edge_key(b, e),
                edge_key(c, e),
                edge_key(d, e),
            }
            missing = specs_missing_edges(adjacency_lists, desired)
            if not missing or len(missing) > 4:
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


def tri_root_hub_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    _seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=6, params=params)
    subset = [int(node) for node in candidates[: min(len(candidates), 8)]]
    specs: list[tuple[float, dict[str, Any]]] = []
    for roots in itertools.combinations(subset[: min(len(subset), 6)], 3):
        remaining = [node for node in subset if node not in roots]
        for center, extra_leaf in itertools.permutations(remaining[: min(len(remaining), 4)], 2):
            a, b, c = [int(value) for value in roots]
            d = int(center)
            e = int(extra_leaf)
            desired = {
                edge_key(d, a),
                edge_key(d, b),
                edge_key(d, c),
                edge_key(d, e),
                edge_key(a, b),
            }
            missing = specs_missing_edges(adjacency_lists, desired)
            if not missing or len(missing) > 4:
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


def three_layer_stacked_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    _seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=5, params=params)
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
        score = sum(float(avg_weights[idx]) for idx in motif_nodes) - 0.07 * len(missing)
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


def hybrid_cross_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    _seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=6, params=params)
    subset = [int(node) for node in candidates[: min(len(candidates), 8)]]
    specs: list[tuple[float, dict[str, Any]]] = []
    for roots in itertools.combinations(subset[: min(len(subset), 6)], 3):
        remaining = [node for node in subset if node not in roots]
        for branches in itertools.combinations(remaining[: min(len(remaining), 5)], 2):
            a, b, c = [int(value) for value in roots]
            d, e = [int(value) for value in branches]
            desired = {
                edge_key(a, b),
                edge_key(b, c),
                edge_key(a, d),
                edge_key(b, d),
                edge_key(b, e),
                edge_key(c, e),
            }
            missing = specs_missing_edges(adjacency_lists, desired)
            if not missing or len(missing) > 5:
                continue
            motif_nodes = sorted({a, b, c, d, e})
            score = sum(float(avg_weights[idx]) for idx in motif_nodes) - 0.07 * len(missing)
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


def build_tri_root_motif_injection(
    adjacency_lists: list[list[int]],
    motif_id: str,
    avg_weights: np.ndarray,
    params: dict[str, Any],
) -> dict[str, Any]:
    if motif_id == 'tri_root_triangle_cluster':
        specs = tri_root_triangle_specs(adjacency_lists, avg_weights, params)
    elif motif_id == 'tri_root_diamond':
        specs = tri_root_diamond_specs(adjacency_lists, avg_weights, params)
    elif motif_id == 'tri_root_hub_micro_star':
        specs = tri_root_hub_specs(adjacency_lists, avg_weights, params)
    elif motif_id == 'three_layer_stacked_motif':
        specs = three_layer_stacked_specs(adjacency_lists, avg_weights, params)
    elif motif_id == 'hybrid_bi_root_tri_root_cross_motif':
        specs = hybrid_cross_specs(adjacency_lists, avg_weights, params)
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
    payload = dict(build_tri_root_motif_injection(adjacency_lists, motif_id, avg_weights, params))
    touched = set(payload_touched_nodes(payload, 'motif'))
    touched.update(int(node) for node in (extra_touched_nodes or []))
    payload['forcing_kind'] = forcing_kind
    payload['forcing_item_id'] = str(forcing_item_id)
    payload['touched_nodes'] = sorted(touched)
    payload['shells'] = compute_shells_from_adjacency(payload['adjacency_lists'], max_shell=int(params.get('max_shell', 4)))
    return payload


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
    if forcing_protocol_id == 'noise_seeded_tri_root_inconsistency':
        return [noise_payloads[noise_id]]
    if forcing_protocol_id == 'tri_root_motif_seeded_inconsistency':
        return [motif_payloads[motif_id]]
    if forcing_protocol_id == 'hybrid_seeded_tri_root_inconsistency':
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
    tri_root_failure_count = int(failure_count if str(payload['forcing_kind']) in {'motif', 'hybrid'} else 0)

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
                'tri_root_motif_correlation': float(motif_corr),
                'kernel_update_consistency_failure_count': int(failure_count),
                'kernel_update_consistency_failure_count_under_tri_root': int(tri_root_failure_count),
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
                'tri_root_forcing_active': int(str(payload['forcing_kind']) in {'motif', 'hybrid'}),
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
        'tri_root_motif_correlation': float(motif_corr),
        'kernel_update_consistency_failure_count': int(template['kernel_update_consistency_failure_count']),
        'kernel_update_consistency_failure_count_under_tri_root': int(template['kernel_update_consistency_failure_count_under_tri_root']),
        'forced_state': template['forced_state'],
        'touched_nodes': list(template['touched_nodes']),
        'motif_nodes': list(template.get('motif_nodes', [])),
        'graph_info': current_graph_info,
        'current_kernel_match': 1,
        'tri_root_forcing_active': int(template['tri_root_forcing_active']),
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
                -int(item['kernel_update_consistency_failure_count_under_tri_root']),
                -int(item['kernel_update_consistency_failure_count']),
                -int(item['candidate']['blocked_merges'] > 0),
                -float(item['mismatch_density']),
                -abs(float(item['tri_root_motif_correlation'])),
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
                -abs(float(item['tri_root_motif_correlation'])),
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
        int(summary['kernel_update_consistency_failure_count_under_tri_root']) > 0
        and int(summary['protected_arrow_step_under_tri_root_inconsistency']) >= 0
        and int(summary['persistent_step']) >= 0
    ):
        return 'tri_root_inconsistency_ordered_persistent'
    if int(summary['kernel_update_consistency_failure_count_under_tri_root']) > 0 and int(summary['protected_arrow_step_under_tri_root_inconsistency']) >= 0:
        return 'tri_root_inconsistency_ordered_transient'
    if (
        int(summary['kernel_update_consistency_failure_count']) > 0
        and int(summary['protected_arrow_step_under_weaker_inconsistency']) >= 0
        and int(summary['persistent_step']) >= 0
    ):
        return 'weaker_inconsistency_ordered_without_tri_root'
    if int(summary['kernel_update_consistency_failure_count']) > 0 and int(summary['protected_arrow_step_under_weaker_inconsistency']) >= 0:
        return 'weaker_inconsistency_transient_without_tri_root'
    if int(summary['kernel_update_consistency_failure_count']) > 0:
        return 'weaker_inconsistency_stable_without_arrow'
    if float(summary['positive_mismatch_density']) > 0.0 and int(summary['protected_arrow_step']) >= 0:
        return 'ordered_without_tri_root_inconsistency'
    if float(summary['positive_mismatch_density']) > 0.0:
        return 'rigid_positive_mismatch'
    return 'stable_without_ordering'


def tri_root_trace_plot(path: Path, run_id: str, trace: list[dict[str, Any]]) -> None:
    steps = [int(item['step_index']) for item in trace]
    mismatch_density = [float(item['mismatch_density']) for item in trace]
    failure_count = [int(item['kernel_update_consistency_failure_count']) for item in trace]
    tri_root_failure_count = [int(item['kernel_update_consistency_failure_count_under_tri_root']) for item in trace]
    blocked = [int(item['protected_blocked_merges']) for item in trace]
    motif_corr = [float(item['tri_root_motif_correlation']) for item in trace]

    fig, axes = plt.subplots(1, 4, figsize=(17.2, 4.2))
    axes[0].plot(steps, [int(item['active_persistence_class_count']) for item in trace], marker='o', color='tab:blue')
    axes[0].set_title('Active persistence classes')
    axes[0].set_xlabel('forced step')
    axes[0].grid(alpha=0.25)

    axes[1].plot(steps, mismatch_density, marker='o', color='tab:red', label='mismatch density')
    axes[1].plot(steps, failure_count, marker='s', color='tab:purple', label='consistency failures')
    axes[1].plot(steps, tri_root_failure_count, marker='^', color='tab:brown', label='tri-root failures')
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
    axes[3].set_title('Tri-root motif correlation')
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

    current_partition = survivor_partition(evaluate_candidate_c08(initial_state, tuple(), seed_edges, 'bounded_path', False, weights, baseline_adj_sets, baseline_shells, seed_shell_distance, params)['final_partition'])
    if not current_partition:
        current_partition = survivor_partition(initial_state['bounded_path_partition'])
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
    arrow_under_weaker: int | None = None
    arrow_under_tri_root: int | None = None
    mismatch_activation_step: int | None = None
    positive_mismatch_steps = 0
    positive_mismatch_density = 0.0
    max_mismatch_density = 0.0
    failure_total = 0
    tri_root_failure_total = 0
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
        tri_root_failure_count = int(chosen['kernel_update_consistency_failure_count_under_tri_root']) if inconsistency_active else 0
        failure_total += failure_count
        tri_root_failure_total += tri_root_failure_count
        failure_max = max(failure_max, failure_count)

        blocked_merges = int(chosen['candidate']['blocked_merges'])
        if blocked_merges > 0:
            arrow_steps.append(int(step_idx + 1))
            if inconsistency_active and failure_count > 0 and arrow_under_weaker is None:
                arrow_under_weaker = int(step_idx + 1)
            if inconsistency_active and tri_root_failure_count > 0 and arrow_under_tri_root is None:
                arrow_under_tri_root = int(step_idx + 1)
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
                'tri_root_motif_correlation': float(chosen['tri_root_motif_correlation']),
                'kernel_update_consistency_failure_count': failure_count,
                'kernel_update_consistency_failure_count_under_tri_root': tri_root_failure_count,
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
        'protected_arrow_step_under_weaker_inconsistency': int(arrow_under_weaker) if arrow_under_weaker is not None else -1,
        'protected_arrow_step_under_tri_root_inconsistency': int(arrow_under_tri_root) if arrow_under_tri_root is not None else -1,
        'persistent_step': int(persistent_step) if persistent_step is not None else -1,
        'final_topology_label': str(current_label),
        'final_active_persistence_class_count': int(trace[-1]['active_persistence_class_count']) if trace else 0,
        'final_seed_component_count': final_seed_components,
        'final_survivor_count': int(len(current_support)),
        'obstruction_rank_final': int(current_state['cycle_rank']),
        'kernel_update_consistency_failure_count': int(failure_total),
        'kernel_update_consistency_failure_count_under_tri_root': int(tri_root_failure_total),
        'kernel_update_consistency_failure_max': int(failure_max),
        'kernel_update_consistency_failure_initial': int(trace[0]['kernel_update_consistency_failure_count']) if trace else 0,
        'kernel_update_consistency_failure_final': int(trace[-1]['kernel_update_consistency_failure_count']) if trace else 0,
        'positive_mismatch_density': float(positive_mismatch_density),
        'max_mismatch_density': float(max_mismatch_density),
        'final_mismatch_density': float(trace[-1]['mismatch_density']) if trace else 0.0,
        'positive_mismatch_steps': int(positive_mismatch_steps),
        'tri_root_motif_correlation': float(max((abs(float(item['tri_root_motif_correlation'])) for item in trace), default=0.0)),
        'irreversibility_scaling_factor': 0.0,
        'protected_blocked_merge_total': int(blocked_total),
        'static_collapse_flag': int(positive_mismatch_density <= mismatch_threshold and (arrow_under_weaker is None)),
        'sequencing_class': '',
        'notes': str(run['notes']),
    }
    summary['sequencing_class'] = sequencing_class(summary)

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_tri_root_trace.png"
    tri_root_trace_plot(trace_path, str(run['run_id']), trace)
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
        'ordered_without_tri_root_inconsistency': 3,
        'weaker_inconsistency_stable_without_arrow': 4,
        'weaker_inconsistency_transient_without_tri_root': 5,
        'weaker_inconsistency_ordered_without_tri_root': 6,
        'tri_root_inconsistency_ordered_transient': 7,
        'tri_root_inconsistency_ordered_persistent': 8,
    }

    matrix_path = WORK_PLOT_DIR / 'stage_c0_9_tri_root_class_matrix.png'
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
    ax.set_title('C0.9 tri-root class matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(matrix_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(matrix_path)

    inconsistency_path = WORK_PLOT_DIR / 'stage_c0_9_tri_root_inconsistency_panel.png'
    x = np.arange(len(reps))
    width = 0.22
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 4.6))
    for idx, protocol in enumerate(protocols):
        fail_values = [int(row_map[(rep, protocol)]['kernel_update_consistency_failure_count']) for rep in reps]
        tri_fail_values = [int(row_map[(rep, protocol)]['kernel_update_consistency_failure_count_under_tri_root']) for rep in reps]
        fail_max = [int(row_map[(rep, protocol)]['kernel_update_consistency_failure_max']) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, fail_values, width=width, label=protocol)
        axes[1].bar(x + (idx - 1) * width, tri_fail_values, width=width, label=protocol)
        axes[2].bar(x + (idx - 1) * width, fail_max, width=width, label=protocol)
    axes[0].set_title('Consistency failure total')
    axes[1].set_title('Tri-root failure total')
    axes[2].set_title('Consistency failure max')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
    axes[0].legend(fontsize=8)
    fig.savefig(inconsistency_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(inconsistency_path)

    arrow_path = WORK_PLOT_DIR / 'stage_c0_9_tri_root_arrow_and_motif_panel.png'
    fig, axes = plt.subplots(1, 4, figsize=(19.2, 4.6))
    for idx, protocol in enumerate(protocols):
        arrow_values = [max(int(row_map[(rep, protocol)]['protected_arrow_step_under_weaker_inconsistency']), 0) for rep in reps]
        tri_arrow_values = [max(int(row_map[(rep, protocol)]['protected_arrow_step_under_tri_root_inconsistency']), 0) for rep in reps]
        mismatch_values = [float(row_map[(rep, protocol)]['positive_mismatch_density']) for rep in reps]
        motif_values = [float(row_map[(rep, protocol)]['tri_root_motif_correlation']) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, arrow_values, width=width, label=protocol)
        axes[1].bar(x + (idx - 1) * width, tri_arrow_values, width=width, label=protocol)
        axes[2].bar(x + (idx - 1) * width, mismatch_values, width=width, label=protocol)
        axes[3].bar(x + (idx - 1) * width, motif_values, width=width, label=protocol)
    axes[0].set_title('Arrow step under weaker inconsistency')
    axes[1].set_title('Arrow step under tri-root inconsistency')
    axes[2].set_title('Positive mismatch density')
    axes[3].set_title('Tri-root motif correlation')
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
    tri_root_failure_runs = sum(1 for row in rows if int(row['kernel_update_consistency_failure_count_under_tri_root']) > 0)
    arrow_runs = sum(1 for row in rows if int(row['protected_arrow_step_under_weaker_inconsistency']) >= 0)
    tri_root_arrow_runs = sum(1 for row in rows if int(row['protected_arrow_step_under_tri_root_inconsistency']) >= 0)
    scaling = float(rows[0]['irreversibility_scaling_factor']) if rows else 0.0
    lines = [
        '# Stage C0.9 Tri-Root Motif Injection v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Architecture notice: {runsheet["architecture_notice"]}',
        f'Baseline notice: {runsheet["baseline_notice"]}',
        '',
        f'Stage summary: explicit consistency failures in `{positive_failure_runs}/{len(rows)}` runs; tri-root-specific failure forcing in `{tri_root_failure_runs}/{len(rows)}` runs; protected arrow under weaker inconsistency in `{arrow_runs}/{len(rows)}` runs; protected arrow under tri-root inconsistency in `{tri_root_arrow_runs}/{len(rows)}` runs; irreversibility scaling factor `{scaling:.2f}` relative to the C0.8 `2/9` baseline.',
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
                f"  - tri-root failure total: `{row['kernel_update_consistency_failure_count_under_tri_root']}`",
                f"  - protected arrow step under weaker inconsistency: `{row['protected_arrow_step_under_weaker_inconsistency']}`",
                f"  - protected arrow step under tri-root inconsistency: `{row['protected_arrow_step_under_tri_root_inconsistency']}`",
                f"  - positive mismatch density: `{float(row['positive_mismatch_density']):.4f}`",
                f"  - tri-root motif correlation: `{float(row['tri_root_motif_correlation']):.4f}`",
            ]
        )
    lines.extend(
        [
            '',
            'Interpretation boundary:',
            '- this scan asks whether tri-root and 3-layer motif forcing can scale the weaker-kernel inconsistency signal opened in C0.8 into a broadly reproducible protected irreversibility mechanism',
            '- topology labels remain persistence-class observables and the protected topology-stability rule stays active',
            '- if consistency failures and protected arrow events both cross the `6/9` bar, the branch supports a broad irreversibility claim; otherwise the correct outcome is a stronger boundary statement, not a forced Phase II interpretation',
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

    tri_root_arrow_runs = sum(1 for row in csv_rows if int(row['protected_arrow_step_under_tri_root_inconsistency']) >= 0)
    scaling = float(tri_root_arrow_runs / 2.0) if csv_rows else 0.0
    for row in csv_rows:
        row['irreversibility_scaling_factor'] = scaling
    for result in results:
        result['summary']['irreversibility_scaling_factor'] = scaling

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
        experiment_slug=str(runsheet.get('experiment_slug', 'stage_c0_9_tri_root_motif_injection')),
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
