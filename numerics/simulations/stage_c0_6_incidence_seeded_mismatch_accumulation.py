#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, load_stage10_defaults, plt, save_atlas_payload
from stage_c0_1_combinatorial_kernel_probe import build_line_graph_cache
from stage_c0_3_incidence_noise_robustness import build_incidence_perturbation, compute_shells_from_adjacency
from stage_c0_4_controlled_motif_injection import build_motif_injection
from stage_c0_5_emergent_compression_preorder import (
    adjacency_sets,
    analyze_support,
    evaluate_candidate,
    initial_edge_weights,
    multi_source_shell_distance,
    next_partition_for_kernel,
    select_active_support,
    top_seed_edges,
    survivor_partition,
    survivor_support,
)
from stage_c0_5_emergent_sequencing import (
    count_seed_components,
    load_reference_rows,
    sequencing_class,
    topology_label_for_support,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_6_incidence_seeded_mismatch_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_6_incidence_seeded_mismatch'
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
    'triangle_cluster',
    'diamond',
    'hub_micro_star',
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
    'initial_active_support',
    'initial_seed_count',
    'steps_completed',
    'tau_sequencing_depth',
    'mismatch_quiet_step',
    'protected_arrow_step',
    'protected_arrow_step_distribution',
    'persistent_step',
    'final_topology_label',
    'final_active_persistence_class_count',
    'final_seed_component_count',
    'final_survivor_count',
    'obstruction_rank_final',
    'kernel_consistency_failures_initial',
    'kernel_consistency_failures_final',
    'kernel_consistency_failures_max',
    'max_mismatch_density',
    'positive_mismatch_density',
    'final_mismatch_density',
    'positive_mismatch_steps',
    'protected_blocked_merge_total',
    'static_collapse_flag',
    'sequencing_class',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.6 incidence-seeded mismatch scan.')
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
        'noise_seeded_protected': 0,
        'motif_seeded_protected': 1,
        'hybrid_seeded_protected': 2,
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


def build_forcing_payloads(
    data: Any,
    representative_id: str,
    weights: np.ndarray,
    params: dict[str, Any],
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    max_shell = int(params.get('max_shell', 3))
    noise_payloads: dict[str, dict[str, Any]] = {}
    for perturbation_id in NOISE_IDS:
        payload = dict(build_incidence_perturbation(data, perturbation_id, params))
        payload['forcing_kind'] = 'noise'
        payload['forcing_item_id'] = str(perturbation_id)
        payload['touched_nodes'] = payload_touched_nodes(payload, 'noise')
        payload['shells'] = compute_shells_from_adjacency(payload['adjacency_lists'], max_shell=max_shell)
        noise_payloads[str(perturbation_id)] = payload

    avg_weights = np.asarray(weights, dtype=float)
    avg_weights = avg_weights / max(float(np.sum(avg_weights)), 1.0e-12)
    motif_payloads: dict[str, dict[str, Any]] = {}
    for motif_id in MOTIF_IDS:
        payload = dict(build_motif_injection(data, representative_id, motif_id, avg_weights, params))
        payload['forcing_kind'] = 'motif'
        payload['forcing_item_id'] = str(motif_id)
        payload['touched_nodes'] = payload_touched_nodes(payload, 'motif')
        payload['shells'] = compute_shells_from_adjacency(payload['adjacency_lists'], max_shell=max_shell)
        motif_payloads[str(motif_id)] = payload
    return noise_payloads, motif_payloads


def forcing_candidates_for_step(
    forcing_protocol_id: str,
    step_idx: int,
    noise_payloads: dict[str, dict[str, Any]],
    motif_payloads: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    noise_id = NOISE_IDS[min(step_idx, len(NOISE_IDS) - 1)]
    motif_id = MOTIF_IDS[step_idx % len(MOTIF_IDS)]
    if forcing_protocol_id == 'noise_seeded_protected':
        return [noise_payloads[noise_id]]
    if forcing_protocol_id == 'motif_seeded_protected':
        return [motif_payloads[motif_id]]
    if forcing_protocol_id == 'hybrid_seeded_protected':
        return [noise_payloads[noise_id], motif_payloads[motif_id]]
    raise ValueError(f'unsupported forcing protocol id: {forcing_protocol_id}')


def forced_mismatch_ledger(
    prev_state: dict[str, Any],
    prev_support: list[int],
    next_state: dict[str, Any],
    next_support: list[int],
    touched_nodes: list[int],
) -> list[int]:
    ledger: list[int] = []
    union_nodes = sorted(
        set(int(node) for node in prev_support)
        .union(int(node) for node in next_support)
        .union(int(node) for node in touched_nodes)
    )
    prev_restricted = dict(prev_state.get('restricted', {}))
    next_restricted = dict(next_state.get('restricted', {}))
    for node in union_nodes:
        prev_sig = prev_state.get('motif_signatures', {}).get(int(node))
        next_sig = next_state.get('motif_signatures', {}).get(int(node))
        prev_degree = len(prev_restricted.get(int(node), set()))
        next_degree = len(next_restricted.get(int(node), set()))
        if prev_sig != next_sig or prev_degree != next_degree or int(node) in touched_nodes:
            ledger.append(int(node))
    return ledger


def sequencing_depth(trace: list[dict[str, Any]]) -> int:
    counts = [int(item['active_persistence_class_count']) for item in trace]
    for idx in range(len(counts)):
        if len(set(counts[idx:])) == 1:
            return int(trace[idx]['step_index'])
    return int(trace[-1]['step_index']) if trace else 0


def sequencing_trace_plot(path: Path, run_id: str, trace: list[dict[str, Any]]) -> None:
    steps = [int(item['step_index']) for item in trace]
    mismatch_density = [float(item['mismatch_density']) for item in trace]
    positive_density = [float(item['positive_mismatch_density']) for item in trace]
    consistency = [int(item['kernel_consistency_failure']) for item in trace]
    blocked = [int(item['protected_blocked_merges']) for item in trace]

    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.2))
    axes[0].plot(steps, [int(item['active_persistence_class_count']) for item in trace], marker='o', color='tab:blue')
    axes[0].set_title('Active persistence classes')
    axes[0].set_xlabel('forced step')
    axes[0].grid(alpha=0.25)

    axes[1].plot(steps, mismatch_density, marker='o', color='tab:red', label='mismatch density')
    axes[1].plot(steps, positive_density, marker='s', color='tab:orange', label='positive mismatch density')
    axes[1].set_title('Mismatch forcing')
    axes[1].set_xlabel('forced step')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)

    axes[2].plot(steps, consistency, marker='o', color='tab:purple', label='kernel consistency failure')
    axes[2].bar(steps, blocked, alpha=0.35, color='tab:gray', label='protected blocks')
    axes[2].set_title('Protected arrow forcing')
    axes[2].set_xlabel('forced step')
    axes[2].grid(alpha=0.25)
    axes[2].legend(fontsize=8)

    fig.suptitle(run_id)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def survivor_overlay_plot(
    path: Path,
    run_id: str,
    midpoints: np.ndarray,
    weights: np.ndarray,
    initial_support: list[int],
    final_support: list[int],
    seed_edges: list[int],
    boundary_type: str,
) -> None:
    coords = np.asarray(midpoints, dtype=float)
    distances = np.abs(coords[:, 2] - 0.5)
    threshold = np.min(distances) + 1.0e-9
    mask = distances <= threshold

    initial_mask = np.zeros(len(coords), dtype=bool)
    initial_mask[np.asarray(initial_support, dtype=int)] = True
    final_mask = np.zeros(len(coords), dtype=bool)
    if final_support:
        final_mask[np.asarray(final_support, dtype=int)] = True
    seed_mask = np.zeros(len(coords), dtype=bool)
    if seed_edges:
        seed_mask[np.asarray(seed_edges, dtype=int)] = True

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    sc = ax.scatter(coords[mask, 0], coords[mask, 1], c=weights[mask], s=14, cmap='Greys', alpha=0.25)
    if np.any(mask & initial_mask):
        ax.scatter(coords[mask & initial_mask, 0], coords[mask & initial_mask, 1], c=weights[mask & initial_mask], s=24, cmap='viridis', alpha=0.65, label='initial support')
    if np.any(mask & final_mask):
        ax.scatter(coords[mask & final_mask, 0], coords[mask & final_mask, 1], s=52, facecolors='none', edgecolors='tab:red', linewidths=1.3, label='final support')
    if np.any(mask & seed_mask):
        ax.scatter(coords[mask & seed_mask, 0], coords[mask & seed_mask, 1], marker='x', s=70, c='tab:orange', linewidths=1.8, label='seed edges')
    ax.set_title(f'Support overlay: {run_id}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend(fontsize=8)
    plt.colorbar(sc, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def evaluate_forcing_option(
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
    forced_state = analyze_support(current_support, weights, full_adj_sets, full_shells, seed_shell_distance, params)
    focus_nodes = [int(node) for node in seed_edges if int(node) in current_support] or list(current_support)

    adjacency_candidate = evaluate_candidate(
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
    shell_candidate = evaluate_candidate(
        forced_state,
        current_partition,
        focus_nodes,
        'graph_shell',
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
    shell_label = topology_label_for_support(
        shell_candidate['next_support'],
        shell_candidate['final_partition'],
        shell_candidate['next_state'],
        seed_edges,
        full_adj_sets,
        protected_label,
    )
    consistency_failure = int(adjacency_label != shell_label)

    results: list[dict[str, Any]] = []
    for kernel_name, candidate, label in (
        ('adjacency', adjacency_candidate, adjacency_label),
        ('graph_shell', shell_candidate, shell_label),
    ):
        ledger = forced_mismatch_ledger(
            prev_state=current_state,
            prev_support=current_support,
            next_state=candidate['next_state'],
            next_support=candidate['next_support'],
            touched_nodes=list(payload['touched_nodes']),
        )
        density = float(len(ledger) / max(len(adjacency_lists), 1))
        results.append(
            {
                'forcing_kind': str(payload['forcing_kind']),
                'forcing_item_id': str(payload['forcing_item_id']),
                'current_topology_label': str(current_label),
                'adjacency_topology_label': str(adjacency_label),
                'graph_shell_topology_label': str(shell_label),
                'active_persistence_class_count': int(len({current_label, adjacency_label, shell_label})),
                'kernel_name': str(kernel_name),
                'candidate': candidate,
                'chosen_topology_label': str(label),
                'preserves_protected_label': int(label == protected_label),
                'mismatch_ledger': ledger,
                'mismatch_density': density,
                'consistency_failure': consistency_failure,
                'forced_state': forced_state,
                'touched_nodes': list(payload['touched_nodes']),
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
    return {
        'forcing_kind': str(template['forcing_kind']),
        'forcing_item_id': str(template['forcing_item_id']),
        'current_topology_label': str(current_label),
        'adjacency_topology_label': str(template['adjacency_topology_label']),
        'graph_shell_topology_label': str(template['graph_shell_topology_label']),
        'active_persistence_class_count': int(template['active_persistence_class_count']),
        'kernel_name': str(current_kernel),
        'candidate': {
            'blocked_merges': max(1, int(template['candidate']['blocked_merges'])),
            'next_support': list(current_support),
            'next_state': current_state,
            'final_partition': current_partition,
        },
        'chosen_topology_label': str(current_label),
        'preserves_protected_label': 1,
        'mismatch_ledger': ledger,
        'mismatch_density': density,
        'consistency_failure': int(template['consistency_failure']),
        'forced_state': template['forced_state'],
        'touched_nodes': list(template['touched_nodes']),
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
) -> dict[str, Any]:
    preserving = [item for item in candidates if int(item['preserves_protected_label']) == 1]
    if not preserving:
        collapsed: list[dict[str, Any]] = []
        seen: set[tuple[str, str]] = set()
        for item in candidates:
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
        candidates = collapsed
    else:
        candidates = preserving

    ranked = sorted(
        candidates,
        key=lambda item: (
            -int(item['candidate']['blocked_merges'] > 0),
            -float(item['mismatch_density']),
            -int(item['consistency_failure']),
            int(item['candidate']['next_state']['support_size']),
            -int(item['current_kernel_match']),
        ),
    )
    return ranked[0]


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
    baseline_graph = build_line_graph_cache(data, max_shell=int(params.get('max_shell', 3)))
    baseline_adj_lists = baseline_graph['adjacency_lists']
    baseline_adj_sets = adjacency_sets(baseline_adj_lists)
    baseline_shells = baseline_graph['shells']
    noise_payloads, motif_payloads = build_forcing_payloads(data, str(run['representative_id']), weights, params)

    initial_support = select_active_support(weights, float(params.get('active_threshold_sigma', 1.0)), int(params.get('min_active_count', 48)))
    seed_edges = top_seed_edges(weights, initial_support, int(params.get('seed_edge_count', 12)))
    seed_shell_distance = multi_source_shell_distance(baseline_adj_lists, seed_edges, int(params.get('seed_shell_cap', 4)))
    initial_state = analyze_support(initial_support, weights, baseline_adj_sets, baseline_shells, seed_shell_distance, params)

    current_partition = survivor_partition(next_partition_for_kernel(initial_state, 'graph_shell'))
    current_support = survivor_support(current_partition)
    current_state = analyze_support(current_support, weights, baseline_adj_sets, baseline_shells, seed_shell_distance, params)
    current_kernel = 'graph_shell'
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
    mismatch_quiet_step: int | None = None
    positive_seen = False
    persistent_step: int | None = None
    last_cycle_rank: int | None = None
    stable_cycle_count = 0
    blocked_total = 0
    positive_mismatch_steps = 0
    max_mismatch_density = 0.0
    positive_mismatch_density = 0.0

    for step_idx in range(max_steps):
        options = forcing_candidates_for_step(str(run['forcing_protocol_id']), step_idx, noise_payloads, motif_payloads)
        candidate_records: list[dict[str, Any]] = []
        for payload in options:
            candidate_records.extend(
                evaluate_forcing_option(
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
        )
        mismatch_density = float(chosen['mismatch_density'])
        max_mismatch_density = max(max_mismatch_density, mismatch_density)
        if mismatch_density > mismatch_threshold:
            positive_seen = True
            positive_mismatch_steps += 1
            positive_mismatch_density = max(positive_mismatch_density, mismatch_density)
        elif positive_seen and mismatch_quiet_step is None:
            mismatch_quiet_step = int(step_idx)

        blocked_merges = int(chosen['candidate']['blocked_merges'])
        if blocked_merges > 0:
            arrow_steps.append(int(step_idx + 1))
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
                'graph_shell_topology_label': str(chosen['graph_shell_topology_label']),
                'active_persistence_class_count': int(chosen['active_persistence_class_count']),
                'chosen_kernel': str(chosen['kernel_name']),
                'chosen_topology_label': str(chosen['chosen_topology_label']),
                'mismatch_ledger_size': int(len(chosen['mismatch_ledger'])),
                'mismatch_density': mismatch_density,
                'positive_mismatch_density': mismatch_density if mismatch_density > mismatch_threshold else 0.0,
                'kernel_consistency_failure': int(chosen['consistency_failure']),
                'protected_blocked_merges': blocked_merges,
                'obstruction_rank_after': cycle_rank,
                'survivor_count_after': int(len(chosen['candidate']['next_support'])),
                'seed_component_count_after': int(count_seed_components(chosen['candidate']['next_support'], seed_edges, chosen['graph_info']['full_adj_sets'])),
            }
        )

        current_support = list(chosen['candidate']['next_support'])
        current_partition = tuple(tuple(int(node) for node in cls) for cls in chosen['candidate']['final_partition'])
        current_state = chosen['candidate']['next_state']
        current_label = str(chosen['chosen_topology_label'])
        current_kernel = str(chosen['kernel_name'])
        current_adj_sets = chosen['graph_info']['full_adj_sets']
        current_graph_info = chosen['graph_info']

    if mismatch_quiet_step is None:
        mismatch_quiet_step = 0 if not positive_seen else -1

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
        'initial_active_support': int(len(initial_support)),
        'initial_seed_count': int(len(seed_edges)),
        'steps_completed': int(len(trace)),
        'tau_sequencing_depth': int(tau_depth),
        'mismatch_quiet_step': int(mismatch_quiet_step),
        'protected_arrow_step': int(arrow_steps[0]) if arrow_steps else -1,
        'protected_arrow_step_distribution': json.dumps(arrow_steps),
        'persistent_step': int(persistent_step) if persistent_step is not None else -1,
        'final_topology_label': str(current_label),
        'final_active_persistence_class_count': int(trace[-1]['active_persistence_class_count']) if trace else 0,
        'final_seed_component_count': final_seed_components,
        'final_survivor_count': int(len(current_support)),
        'obstruction_rank_final': int(current_state['cycle_rank']),
        'kernel_consistency_failures_initial': int(trace[0]['kernel_consistency_failure']) if trace else 0,
        'kernel_consistency_failures_final': int(trace[-1]['kernel_consistency_failure']) if trace else 0,
        'kernel_consistency_failures_max': int(max((item['kernel_consistency_failure'] for item in trace), default=0)),
        'max_mismatch_density': float(max_mismatch_density),
        'positive_mismatch_density': float(positive_mismatch_density),
        'final_mismatch_density': float(trace[-1]['mismatch_density']) if trace else 0.0,
        'positive_mismatch_steps': int(positive_mismatch_steps),
        'protected_blocked_merge_total': int(blocked_total),
        'static_collapse_flag': int(positive_mismatch_density <= mismatch_threshold and not arrow_steps),
        'sequencing_class': '',
        'notes': str(run['notes']),
    }
    summary['sequencing_class'] = sequencing_class(summary)

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_forced_trace.png"
    sequencing_trace_plot(trace_path, str(run['run_id']), trace)
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
        'protected_ordered_transient': 2,
        'protected_ordered_persistent': 3,
    }

    matrix_path = WORK_PLOT_DIR / 'stage_c0_6_sequencing_class_matrix.png'
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
    ax.set_title('C0.6 sequencing class matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(matrix_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(matrix_path)

    mismatch_path = WORK_PLOT_DIR / 'stage_c0_6_positive_mismatch_panel.png'
    x = np.arange(len(reps))
    width = 0.22
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6))
    for idx, protocol in enumerate(protocols):
        density_values = [float(row_map[(rep, protocol)]['positive_mismatch_density']) for rep in reps]
        step_values = [int(row_map[(rep, protocol)]['positive_mismatch_steps']) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, density_values, width=width, label=protocol)
        axes[1].bar(x + (idx - 1) * width, step_values, width=width, label=protocol)
    axes[0].set_title('Positive mismatch density')
    axes[1].set_title('Positive mismatch steps')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
    axes[0].legend(fontsize=8)
    fig.savefig(mismatch_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(mismatch_path)

    arrow_path = WORK_PLOT_DIR / 'stage_c0_6_protected_arrow_panel.png'
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6))
    for idx, protocol in enumerate(protocols):
        arrow_values = [max(int(row_map[(rep, protocol)]['protected_arrow_step']), 0) for rep in reps]
        blocked_values = [int(row_map[(rep, protocol)]['protected_blocked_merge_total']) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, arrow_values, width=width, label=protocol)
        axes[1].bar(x + (idx - 1) * width, blocked_values, width=width, label=protocol)
    axes[0].set_title('First protected arrow step')
    axes[1].set_title('Protected blocked merges')
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
    distribution: Counter[int] = Counter()
    for row in rows:
        for step in json.loads(str(row['protected_arrow_step_distribution'])):
            distribution[int(step)] += 1
    positive_runs = sum(1 for row in rows if float(row['positive_mismatch_density']) > 0.0)
    arrow_runs = sum(1 for row in rows if int(row['protected_arrow_step']) >= 0)

    lines = [
        '# Stage C0.6 Incidence-Seeded Mismatch Accumulation v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Architecture notice: {runsheet["architecture_notice"]}',
        f'Baseline notice: {runsheet["baseline_notice"]}',
        '',
        f'Stage summary: positive mismatch in `{positive_runs}/{len(rows)}` runs; protected arrow forcing in `{arrow_runs}/{len(rows)}` runs.',
        f'Protected arrow step distribution: `{dict(sorted(distribution.items()))}`',
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
                f"  - positive mismatch density: `{float(row['positive_mismatch_density']):.4f}`",
                f"  - positive mismatch steps: `{row['positive_mismatch_steps']}`",
                f"  - protected arrow step distribution: `{row['protected_arrow_step_distribution']}`",
                f"  - kernel consistency max: `{row['kernel_consistency_failures_max']}`",
            ]
        )
    lines.extend(
        [
            '',
            'Interpretation boundary:',
            '- this scan asks whether external combinatorial forcing can produce positive mismatch density while the protected topology-stability rule remains active',
            '- topology labels are still treated as persistence-class observables on the no-distance branch',
            '- if positive mismatch remains absent or isolated, the branch is still a rigid persistence-class control rather than an incidence-driven sequencing mechanism',
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
        experiment_slug=str(runsheet.get('experiment_slug', 'stage_c0_6_incidence_seeded_mismatch_accumulation')),
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
