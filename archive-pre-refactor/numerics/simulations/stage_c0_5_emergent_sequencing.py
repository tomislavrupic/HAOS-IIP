#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, load_stage10_defaults, plt, save_atlas_payload
from stage_c0_1_combinatorial_kernel_probe import build_line_graph_cache
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

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_5_emergent_sequencing_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_5_emergent_sequencing'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'representative_id',
    'sequencing_mode_id',
    'sequencing_mode_label',
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
    'final_mismatch_density',
    'protected_blocked_merge_total',
    'static_collapse_flag',
    'sequencing_class',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.5 emergent sequencing scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_runs(runs: list[dict[str, Any]], run_ids: list[str]) -> list[dict[str, Any]]:
    if not run_ids:
        return runs
    wanted = set(run_ids)
    return [run for run in runs if run['run_id'] in wanted]


def lookup_by_id(items: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    return {str(item[key]): item for item in items}


def sequencing_mode_order(mode_id: str) -> int:
    return {
        'adjacency_sequence': 0,
        'graph_shell_sequence': 1,
        'protected_alternating_sequence': 2,
    }.get(mode_id, 99)


def latest_csv(slug: str) -> Path | None:
    candidates = sorted((REPO_ROOT / 'data').glob(f'*_{slug}.csv'))
    return candidates[-1] if candidates else None


def load_reference_rows() -> tuple[dict[str, str], dict[tuple[str, str], dict[str, str]]]:
    protected_labels: dict[str, str] = {}
    motif_csv = latest_csv('stage_c0_4_controlled_motif_injection')
    if motif_csv is not None:
        with motif_csv.open('r', encoding='utf-8', newline='') as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rep = str(row['representative_id'])
                protected_labels.setdefault(rep, str(row['baseline_topology']))

    c01_ref: dict[tuple[str, str], dict[str, str]] = {}
    c01_csv = latest_csv('stage_c0_1_combinatorial_kernel_probe')
    if c01_csv is not None:
        with c01_csv.open('r', encoding='utf-8', newline='') as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                key = (str(row['representative_id']), str(row['kernel_class_id']))
                c01_ref[key] = {
                    'interaction_label': str(row['interaction_label']),
                    'coarse_label': str(row['coarse_label']),
                }
    return protected_labels, c01_ref


def count_seed_components(support: list[int], seed_edges: list[int], full_adj_sets: list[set[int]]) -> int:
    support_set = set(int(node) for node in support)
    seed_set = {int(node) for node in seed_edges if int(node) in support_set}
    if not seed_set:
        return 0
    visited: set[int] = set()
    components = 0
    for seed in sorted(seed_set):
        if seed in visited:
            continue
        components += 1
        stack = [seed]
        visited.add(seed)
        while stack:
            cur = stack.pop()
            for nbr in full_adj_sets[cur]:
                if nbr in seed_set and nbr in support_set and nbr not in visited:
                    visited.add(nbr)
                    stack.append(nbr)
    return components


def topology_label_for_support(
    support: list[int],
    partition: tuple[tuple[int, ...], ...],
    state: dict[str, Any],
    seed_edges: list[int],
    full_adj_sets: list[set[int]],
    fallback_label: str,
) -> str:
    if not support:
        return fallback_label
    seed_components = count_seed_components(support, seed_edges, full_adj_sets)
    class_sizes = sorted((len(cls) for cls in partition), reverse=True)
    largest_class = max(class_sizes, default=0)
    class_count = len(class_sizes)
    cycle_rank = int(state['cycle_rank'])

    if cycle_rank > 0 and largest_class <= 4:
        return 'trapped_local'
    if seed_components == 3 and cycle_rank == 0 and largest_class >= 6:
        return 'split_channel_exchange'
    if seed_components >= 4 and cycle_rank == 0 and largest_class <= 4 and class_count >= 5:
        return 'trapped_local'
    if seed_components <= 2:
        return 'smeared_transfer'
    if seed_components >= 4 and cycle_rank == 0 and class_count <= 3:
        return 'smeared_transfer'
    return fallback_label


def incidence_mismatch_ledger(
    prev_state: dict[str, Any],
    prev_support: list[int],
    next_state: dict[str, Any],
    next_support: list[int],
    protected_label: str,
    candidate_label: str,
) -> list[int]:
    if candidate_label == protected_label:
        return []
    union_nodes = sorted(set(int(node) for node in prev_support).union(int(node) for node in next_support))
    ledger: list[int] = []
    for node in union_nodes:
        prev_sig = prev_state['motif_signatures'].get(int(node))
        next_sig = next_state['motif_signatures'].get(int(node))
        if prev_sig != next_sig:
            ledger.append(int(node))
    return ledger


def sequencing_class(summary: dict[str, Any]) -> str:
    if int(summary['static_collapse_flag']) == 1 and int(summary['protected_arrow_step']) < 0:
        return 'static_regime'
    if int(summary['protected_arrow_step']) >= 0 and int(summary['persistent_step']) >= 0:
        return 'protected_ordered_persistent'
    if int(summary['protected_arrow_step']) >= 0:
        return 'protected_ordered_transient'
    return 'stable_without_ordering'


def sequencing_trace_plot(path: Path, run_id: str, trace: list[dict[str, Any]]) -> None:
    steps = [int(item['step_index']) for item in trace]
    class_counts = [int(item['active_persistence_class_count']) for item in trace]
    mismatch_density = [float(item['mismatch_density']) for item in trace]
    consistency = [int(item['kernel_consistency_failure']) for item in trace]
    obstruction_rank = [int(item['obstruction_rank_after']) for item in trace]
    blocked = [int(item['protected_blocked_merges']) for item in trace]

    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.3))
    axes[0].plot(steps, class_counts, marker='o', color='tab:blue')
    axes[0].set_title('Active persistence classes')
    axes[0].set_xlabel('sequencing step')
    axes[0].grid(alpha=0.25)

    axes[1].plot(steps, mismatch_density, marker='o', color='tab:red', label='mismatch density')
    axes[1].plot(steps, consistency, marker='s', color='tab:purple', label='kernel consistency failure')
    axes[1].set_title('Mismatch accumulation')
    axes[1].set_xlabel('sequencing step')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)

    axes[2].plot(steps, obstruction_rank, marker='o', color='tab:green', label='obstruction rank')
    axes[2].bar(steps, blocked, alpha=0.35, color='tab:gray', label='protected blocks')
    axes[2].set_title('Protected stability and obstruction rank')
    axes[2].set_xlabel('sequencing step')
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
    slice_value = 0.5 if boundary_type == 'periodic' else 0.5
    distances = np.abs(coords[:, 2] - slice_value)
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

    fig, ax = plt.subplots(figsize=(6.1, 5.1))
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


def sequencing_depth(trace: list[dict[str, Any]]) -> int:
    counts = [int(item['active_persistence_class_count']) for item in trace]
    for idx in range(len(counts)):
        if len(set(counts[idx:])) == 1:
            return int(trace[idx]['step_index'])
    return int(trace[-1]['step_index']) if trace else 0


def simulate_run(
    run: dict[str, Any],
    representative: dict[str, Any],
    mode: dict[str, Any],
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
    line_graph = build_line_graph_cache(data, max_shell=int(params.get('max_shell', 3)))
    full_adj_lists = line_graph['adjacency_lists']
    full_adj_sets = adjacency_sets(full_adj_lists)
    full_shells = line_graph['shells']

    initial_support = select_active_support(weights, float(params.get('active_threshold_sigma', 1.0)), int(params.get('min_active_count', 48)))
    seed_edges = top_seed_edges(weights, initial_support, int(params.get('seed_edge_count', 12)))
    seed_shell_distance = multi_source_shell_distance(full_adj_lists, seed_edges, int(params.get('seed_shell_cap', 4)))

    initial_state = analyze_support(initial_support, weights, full_adj_sets, full_shells, seed_shell_distance, params)
    start_kernel = 'adjacency' if str(run['sequencing_mode_id']) == 'adjacency_sequence' else 'graph_shell'
    current_partition = survivor_partition(next_partition_for_kernel(initial_state, start_kernel))
    current_support = survivor_support(current_partition)
    current_state = analyze_support(current_support, weights, full_adj_sets, full_shells, seed_shell_distance, params)
    current_kernel = start_kernel
    current_label = topology_label_for_support(current_support, current_partition, current_state, seed_edges, full_adj_sets, protected_label)
    current_focus_nodes = [int(node) for node in seed_edges if int(node) in current_support] or list(current_support)

    trace: list[dict[str, Any]] = []
    arrow_step: int | None = None
    mismatch_quiet_step: int | None = None
    persistent_step: int | None = None
    blocked_total = 0
    last_cycle_rank: int | None = None
    stable_cycle_count = 0

    mismatch_threshold = float(params.get('mismatch_density_threshold', 0.001))
    rho_window = int(params.get('rho_stability_window', 2))
    max_steps = int(params.get('max_steps', 8))

    for step_idx in range(0, max_steps + 1):
        adjacency_candidate = evaluate_candidate(
            current_state,
            current_partition,
            current_focus_nodes,
            'adjacency',
            str(run['sequencing_mode_id']) == 'protected_alternating_sequence',
            weights,
            full_adj_sets,
            full_shells,
            seed_shell_distance,
            params,
        )
        shell_candidate = evaluate_candidate(
            current_state,
            current_partition,
            current_focus_nodes,
            'graph_shell',
            str(run['sequencing_mode_id']) == 'protected_alternating_sequence',
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
        active_labels = sorted({current_label, adjacency_label, shell_label})
        consistency_failure = int(adjacency_label != shell_label)

        if step_idx == max_steps:
            chosen_candidate = {
                'next_support': current_support,
                'next_state': current_state,
                'final_partition': current_partition,
                'blocked_merges': 0,
                'kernel_name': current_kernel,
            }
            chosen_label = current_label
            ledger = []
        else:
            adjacency_ledger = incidence_mismatch_ledger(current_state, current_support, adjacency_candidate['next_state'], adjacency_candidate['next_support'], protected_label, adjacency_label)
            shell_ledger = incidence_mismatch_ledger(current_state, current_support, shell_candidate['next_state'], shell_candidate['next_support'], protected_label, shell_label)
            adjacency_density = float(len(adjacency_ledger) / max(len(full_adj_lists), 1))
            shell_density = float(len(shell_ledger) / max(len(full_adj_lists), 1))

            if str(run['sequencing_mode_id']) == 'adjacency_sequence':
                chosen_candidate = adjacency_candidate
                chosen_label = adjacency_label
                ledger = adjacency_ledger
            elif str(run['sequencing_mode_id']) == 'graph_shell_sequence':
                chosen_candidate = shell_candidate
                chosen_label = shell_label
                ledger = shell_ledger
            else:
                ranked = sorted(
                    [
                        (adjacency_candidate, adjacency_label, adjacency_ledger, adjacency_density, 'adjacency'),
                        (shell_candidate, shell_label, shell_ledger, shell_density, 'graph_shell'),
                    ],
                    key=lambda item: (
                        float(item[3]),
                        int(item[0]['blocked_merges']),
                        int(item[0]['next_state']['support_size']),
                        0 if item[4] == current_kernel else 1,
                    ),
                )
                chosen_candidate, chosen_label, ledger, _density, _name = ranked[0]
                if int(chosen_candidate['blocked_merges']) > 0 and arrow_step is None:
                    arrow_step = int(step_idx + 1)
                blocked_total += int(chosen_candidate['blocked_merges'])

        mismatch_density = float(len(ledger) / max(len(full_adj_lists), 1))
        if mismatch_quiet_step is None and mismatch_density <= mismatch_threshold:
            mismatch_quiet_step = int(step_idx)

        cycle_rank = int(chosen_candidate['next_state']['cycle_rank'])
        if last_cycle_rank == cycle_rank:
            stable_cycle_count += 1
        else:
            stable_cycle_count = 1
            last_cycle_rank = cycle_rank
        if persistent_step is None and stable_cycle_count >= rho_window:
            persistent_step = int(step_idx)

        trace.append(
            {
                'step_index': int(step_idx),
                'current_topology_label': current_label,
                'adjacency_topology_label': adjacency_label,
                'graph_shell_topology_label': shell_label,
                'active_persistence_class_count': int(len(active_labels)),
                'kernel_consistency_failure': consistency_failure,
                'chosen_kernel': str(chosen_candidate['kernel_name']),
                'chosen_topology_label': chosen_label,
                'mismatch_ledger_size': int(len(ledger)),
                'mismatch_density': mismatch_density,
                'protected_blocked_merges': int(chosen_candidate['blocked_merges']),
                'obstruction_rank_after': cycle_rank,
                'survivor_count_after': int(len(chosen_candidate['next_support'])),
                'seed_component_count_after': int(count_seed_components(chosen_candidate['next_support'], seed_edges, full_adj_sets)),
            }
        )

        if step_idx == max_steps:
            break

        stabilized = (
            chosen_candidate['next_support'] == current_support
            and chosen_candidate['final_partition'] == current_partition
        )
        current_support = list(chosen_candidate['next_support'])
        current_partition = tuple(tuple(int(node) for node in cls) for cls in chosen_candidate['final_partition'])
        current_state = chosen_candidate['next_state']
        current_label = chosen_label
        current_kernel = str(chosen_candidate['kernel_name'])
        current_focus_nodes = [int(node) for node in seed_edges if int(node) in current_support] or list(current_support)
        if stabilized:
            break

    tau_depth = sequencing_depth(trace)
    static_collapse = int(mismatch_quiet_step is not None and (arrow_step is None or int(mismatch_quiet_step) < int(arrow_step)))
    reference_kernel = 'graph_shell' if str(run['sequencing_mode_id']) != 'adjacency_sequence' else 'pure_adjacency'
    reference = c01_ref.get((str(run['representative_id']), reference_kernel), {})
    final_seed_components = int(count_seed_components(current_support, seed_edges, full_adj_sets))
    summary = {
        'run_id': str(run['run_id']),
        'representative_id': str(run['representative_id']),
        'sequencing_mode_id': str(run['sequencing_mode_id']),
        'sequencing_mode_label': str(mode['label']),
        'protected_topology_label': protected_label,
        'reference_interaction_label': str(reference.get('interaction_label', 'unknown')),
        'reference_coarse_label': str(reference.get('coarse_label', 'unknown')),
        'resolution': resolution,
        'initial_active_support': int(len(initial_support)),
        'initial_seed_count': int(len(seed_edges)),
        'steps_completed': int(trace[-1]['step_index']) if trace else 0,
        'tau_sequencing_depth': int(tau_depth),
        'mismatch_quiet_step': int(mismatch_quiet_step) if mismatch_quiet_step is not None else -1,
        'protected_arrow_step': int(arrow_step) if arrow_step is not None else -1,
        'persistent_step': int(persistent_step) if persistent_step is not None else -1,
        'final_topology_label': str(current_label),
        'final_active_persistence_class_count': int(trace[-1]['active_persistence_class_count']) if trace else 0,
        'final_seed_component_count': final_seed_components,
        'final_survivor_count': int(len(current_support)),
        'obstruction_rank_final': int(current_state['cycle_rank']),
        'kernel_consistency_failures_initial': int(trace[0]['kernel_consistency_failure']) if trace else 0,
        'kernel_consistency_failures_final': int(trace[-1]['kernel_consistency_failure']) if trace else 0,
        'kernel_consistency_failures_max': int(max((item['kernel_consistency_failure'] for item in trace), default=0)),
        'max_mismatch_density': float(max((item['mismatch_density'] for item in trace), default=0.0)),
        'final_mismatch_density': float(trace[-1]['mismatch_density']) if trace else 0.0,
        'protected_blocked_merge_total': int(blocked_total),
        'static_collapse_flag': static_collapse,
        'sequencing_class': '',
        'notes': str(run['notes']),
    }
    summary['sequencing_class'] = sequencing_class(summary)

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_sequencing_trace.png"
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
        'sequencing_mode': mode,
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
    modes = sorted({str(row['sequencing_mode_id']) for row in rows}, key=sequencing_mode_order)
    row_map = {(str(row['representative_id']), str(row['sequencing_mode_id'])): row for row in rows}

    class_scores = {
        'static_regime': 0,
        'stable_without_ordering': 1,
        'protected_ordered_transient': 2,
        'protected_ordered_persistent': 3,
    }

    matrix_path = WORK_PLOT_DIR / 'stage_c0_5_sequencing_class_matrix.png'
    matrix = np.zeros((len(reps), len(modes)), dtype=float)
    for i, rep in enumerate(reps):
        for j, mode in enumerate(modes):
            matrix[i, j] = float(class_scores.get(str(row_map[(rep, mode)]['sequencing_class']), 0))
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    im = ax.imshow(matrix, cmap='viridis', vmin=0, vmax=max(class_scores.values()))
    ax.set_xticks(range(len(modes)))
    ax.set_xticklabels(modes, rotation=20, ha='right')
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('Sequencing class matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(matrix_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(matrix_path)

    depth_path = WORK_PLOT_DIR / 'stage_c0_5_sequencing_depth_panel.png'
    x = np.arange(len(reps))
    width = 0.22
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6))
    for idx, mode in enumerate(modes):
        tau_values = [int(row_map[(rep, mode)]['tau_sequencing_depth']) for rep in reps]
        persistent_values = [max(int(row_map[(rep, mode)]['persistent_step']), 0) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, tau_values, width=width, label=mode)
        axes[1].bar(x + (idx - 1) * width, persistent_values, width=width, label=mode)
    axes[0].set_title('Sequencing depth')
    axes[1].set_title('Persistence depth')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
    axes[0].legend(fontsize=8)
    fig.savefig(depth_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(depth_path)

    mismatch_path = WORK_PLOT_DIR / 'stage_c0_5_mismatch_panel.png'
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6))
    for idx, mode in enumerate(modes):
        max_delta = [float(row_map[(rep, mode)]['max_mismatch_density']) for rep in reps]
        arrows = [max(int(row_map[(rep, mode)]['protected_arrow_step']), 0) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, max_delta, width=width, label=mode)
        axes[1].bar(x + (idx - 1) * width, arrows, width=width, label=mode)
    axes[0].set_title('Mismatch density by mode')
    axes[1].set_title('Protected arrow step')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
    axes[0].legend(fontsize=8)
    fig.savefig(mismatch_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(mismatch_path)
    return plot_paths


def write_note(path: Path, json_rel: str, csv_rel: str, rows: list[dict[str, Any]], runsheet: dict[str, Any]) -> None:
    lines = [
        '# Stage C0.5 Emergent Sequencing v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Architecture notice: {runsheet["architecture_notice"]}',
        f'Baseline notice: {runsheet["baseline_notice"]}',
        '',
        'Per-run summary:',
        '',
    ]
    for row in rows:
        lines.extend(
            [
                f"- `{row['sequencing_mode_label']}` / `{row['representative_id']}`",
                f"  - protected topology label: `{row['protected_topology_label']}`",
                f"  - sequencing class: `{row['sequencing_class']}`",
                f"  - sequencing depth: `{row['tau_sequencing_depth']}`",
                f"  - mismatch quiet step: `{row['mismatch_quiet_step']}`",
                f"  - protected arrow step: `{row['protected_arrow_step']}`",
                f"  - persistence step: `{row['persistent_step']}`",
                f"  - final topology label: `{row['final_topology_label']}`",
                f"  - max mismatch density: `{row['max_mismatch_density']:.4f}`",
            ]
        )
    lines.extend(
        [
            '',
            'Interpretation boundary:',
            '- this scan derives sequencing from incidence-mismatch accumulation on combinatorial support updates',
            '- topology labels are treated as persistence-class observables on the no-distance branch',
            '- it does not claim physical time, geometry emergence, particles, or universal irreversibility',
        ]
    )
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    defaults = load_stage10_defaults()
    base = dict(runsheet['base_seed_reference'])
    representatives = lookup_by_id(runsheet['representatives'], 'representative_id')
    modes = lookup_by_id(runsheet['sequencing_modes'], 'sequencing_mode_id')
    protected_labels, c01_ref = load_reference_rows()

    results: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in selected_runs(runsheet['runs'], args.run_ids):
        representative = representatives[str(run['representative_id'])]
        mode = modes[str(run['sequencing_mode_id'])]
        result, run_plots = simulate_run(run, representative, mode, runsheet, defaults, base, protected_labels, c01_ref)
        results.append(result)
        csv_rows.append(dict(result['summary']))
        plot_paths.extend(run_plots)

    csv_rows.sort(key=lambda row: (str(row['representative_id']), sequencing_mode_order(str(row['sequencing_mode_id']))))
    plot_paths.extend(create_summary_plots(csv_rows))

    payload = {
        'stage': str(runsheet['stage']),
        'description': str(runsheet['description']),
        'architecture_notice': str(runsheet['architecture_notice']),
        'baseline_notice': str(runsheet['baseline_notice']),
        'results': results,
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug=str(runsheet.get('experiment_slug', 'stage_c0_5_emergent_sequencing')),
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
