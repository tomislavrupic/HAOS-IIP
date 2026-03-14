#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import (
    ATLAS_NOTES,
    AtlasRun,
    PLOTS,
    REPO_ROOT,
    append_log,
    build_regular_complex,
    build_scalar_initial_state,
    build_transverse_initial_state,
    build_transverse_setup,
    create_field_snapshot,
    create_spectral_trace,
    create_trace_plot,
    load_stage10_defaults,
    plt,
    save_atlas_payload,
    simulate_atlas_run,
    suggested_dt,
)
from stage10_perturbations import build_perturbed_operators
from stage10_regime_labels import classify_run

NOTE_PATH = ATLAS_NOTES / 'Atlas_1_Perturbation_Resilience_v1.md'
RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage10_atlas1_runs.json'
BASELINE_CSV_PATH = REPO_ROOT / 'data' / '20260314_143113_stage10_atlas0_baseline_morphology.csv'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage10_atlas1_plots'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'atlas_phase',
    'baseline_run_id',
    'perturbation_axis',
    'perturbation_strength',
    'operator_sector',
    'boundary_type',
    'random_seed',
    'baseline_regime',
    'perturbed_regime',
    'perturbed_regime_perturbed_constraint',
    'baseline_sector_label',
    'sector_label',
    'sector_label_perturbed_constraint',
    'persistence_score',
    'persistence_score_perturbed_constraint',
    'constraint_max',
    'constraint_max_perturbed',
    'baseline_constraint_max',
    'sector_leakage',
    'baseline_sector_leakage',
    'norm_drift',
    'anisotropy_score',
    'coherence_score',
    'baseline_coherence_score',
    'coherence_delta',
    'center_shift',
    'width_ratio',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 10B Atlas-1 perturbation resilience pilot.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    parser.add_argument('--max-runs', type=int, default=0)
    return parser.parse_args()


def load_baseline_rows(path: Path) -> dict[str, dict[str, Any]]:
    with path.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key in (
            'constraint_max',
            'sector_leakage',
            'norm_drift',
            'anisotropy_score',
            'coherence_score',
            'center_shift',
            'width_ratio',
            'spectral_centroid_final',
            'spectral_spread_final',
            'central_k',
            'bandwidth',
            'amplitude',
        ):
            row[key] = float(row[key])
        row['packet_count'] = int(row['packet_count'])
        row['random_seed'] = int(row['random_seed'])
        row['proto_spacetime_score'] = int(row['proto_spacetime_score'])
    return {row['run_id']: row for row in rows}


def load_runsheet(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding='utf-8'))
    return payload['runs']


def selected_runs(all_runs: list[dict[str, Any]], run_ids: list[str], max_runs: int) -> list[dict[str, Any]]:
    runs = all_runs
    if run_ids:
        wanted = set(run_ids)
        runs = [run for run in runs if run['run_id'] in wanted]
    if max_runs > 0:
        runs = runs[:max_runs]
    return runs


def atlas_run_from_baseline(row: dict[str, Any]) -> AtlasRun:
    return AtlasRun(
        run_id=row['run_id'],
        atlas_phase=row['atlas_phase'],
        graph_type=row['graph_type'],
        kernel_type=row['kernel_type'],
        operator_sector=row['operator_sector'],
        boundary_type=row['boundary_type'],
        initial_seed_type=row['initial_seed_type'],
        central_k=float(row['central_k']),
        bandwidth=float(row['bandwidth']),
        amplitude=float(row['amplitude']),
        phase_pattern=row['phase_pattern'],
        packet_count=int(row['packet_count']),
        random_seed=int(row['random_seed']),
        notes=row['notes'],
    )


def make_leakage_fn(projector):
    def leakage(vec: np.ndarray) -> float:
        denom = max(float(np.linalg.norm(vec)), 1.0e-12)
        return float(np.linalg.norm(np.asarray(projector(vec), dtype=float) - vec) / denom)

    return leakage


def replay_constraint_max(q0: np.ndarray, v0: np.ndarray, dt: float, steps: int, operator, divergence_op, project=None) -> float:
    q = np.asarray(q0, dtype=float).copy()
    v = np.asarray(v0, dtype=float).copy()
    if project is not None:
        q = np.asarray(project(q), dtype=float)
        v = np.asarray(project(v), dtype=float)
    max_constraint = float(np.linalg.norm(divergence_op @ q))
    for _ in range(steps):
        acc = -np.asarray(operator @ q, dtype=float)
        if project is not None:
            acc = np.asarray(project(acc), dtype=float)
        v_half = v + 0.5 * dt * acc
        q = q + dt * v_half
        if project is not None:
            q = np.asarray(project(q), dtype=float)
        acc_new = -np.asarray(operator @ q, dtype=float)
        if project is not None:
            acc_new = np.asarray(project(acc_new), dtype=float)
        v = v_half + 0.5 * dt * acc_new
        if project is not None:
            v = np.asarray(project(v), dtype=float)
        max_constraint = max(max_constraint, float(np.linalg.norm(divergence_op @ q)))
    return max_constraint


def final_values_from_sim(q0: np.ndarray, v0: np.ndarray, dt: float, steps: int, operator, project=None) -> np.ndarray:
    q = np.asarray(q0, dtype=float).copy()
    v = np.asarray(v0, dtype=float).copy()
    if project is not None:
        q = np.asarray(project(q), dtype=float)
        v = np.asarray(project(v), dtype=float)
    for _ in range(steps):
        acc = -np.asarray(operator @ q, dtype=float)
        if project is not None:
            acc = np.asarray(project(acc), dtype=float)
        v_half = v + 0.5 * dt * acc
        q = q + dt * v_half
        if project is not None:
            q = np.asarray(project(q), dtype=float)
        acc_new = -np.asarray(operator @ q, dtype=float)
        if project is not None:
            acc_new = np.asarray(project(acc_new), dtype=float)
        v = v_half + 0.5 * dt * acc_new
        if project is not None:
            v = np.asarray(project(v), dtype=float)
    return q


def create_transition_heatmap(path: Path, rows: list[dict[str, Any]]) -> None:
    order = [
        'Ballistic coherent',
        'Ballistic dispersive',
        'Diffusive',
        'Localized',
        'Oscillatory trapped',
        'Chaotic or irregular',
        'Fragmenting',
        'Metastable structured',
    ]
    index = {label: idx for idx, label in enumerate(order)}
    matrix = np.zeros((len(order), len(order)), dtype=int)
    for row in rows:
        i = index[row['baseline_regime']]
        j = index[row['perturbed_regime']]
        matrix[i, j] += 1
    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    im = ax.imshow(matrix, cmap='Blues')
    ax.set_xticks(range(len(order)), order, rotation=30, ha='right')
    ax.set_yticks(range(len(order)), order)
    ax.set_xlabel('perturbed regime')
    ax.set_ylabel('baseline regime')
    ax.set_title('Atlas-1 transition counts')
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix[i, j]), ha='center', va='center', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_persistence_bar(path: Path, rows: list[dict[str, Any]]) -> None:
    order = ['Ballistic coherent', 'Ballistic dispersive', 'Diffusive', 'Chaotic or irregular', 'Fragmenting']
    values = []
    for label in order:
        sub = [row['persistence_score'] for row in rows if row['baseline_regime'] == label]
        values.append(float(np.mean(sub)) if sub else 0.0)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.bar(range(len(order)), values, color='tab:blue', alpha=0.8)
    ax.set_xticks(range(len(order)), order, rotation=30, ha='right')
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel('mean persistence score')
    ax.set_title('Atlas-1 persistence by baseline regime')
    ax.grid(axis='y', alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_coherence_change_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    colors = {
        'edge_weight_noise': 'tab:blue',
        'sparse_graph_defects': 'tab:orange',
        'anisotropic_bias': 'tab:green',
        'kernel_width_jitter': 'tab:red',
    }
    fig, ax = plt.subplots(figsize=(6.1, 4.8))
    for axis, color in colors.items():
        sub = [row for row in rows if row['perturbation_axis'] == axis]
        ax.scatter(
            [row['baseline_coherence_score'] for row in sub],
            [row['coherence_score'] for row in sub],
            label=axis,
            color=color,
            alpha=0.8,
        )
    lims = [0.35, 1.0]
    ax.plot(lims, lims, color='black', linestyle='--', linewidth=1.0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel('baseline coherence score')
    ax.set_ylabel('perturbed coherence score')
    ax.set_title('Atlas-1 coherence response')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_constraint_leakage_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    axes = ['edge_weight_noise', 'sparse_graph_defects', 'anisotropic_bias', 'kernel_width_jitter']
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6.4, 8.6), sharex=True)
    for idx, axis in enumerate(axes):
        sub = [row for row in rows if row['perturbation_axis'] == axis]
        xs = np.full(len(sub), idx, dtype=float)
        jitter = np.linspace(-0.12, 0.12, num=max(len(sub), 1))
        xs = xs + jitter[:len(sub)]
        ax1.scatter(xs, [max(row['constraint_max'], 1.0e-18) for row in sub], alpha=0.8)
        ax2.scatter(xs, [max(row['constraint_max_perturbed'], 1.0e-18) for row in sub], alpha=0.8)
        ax3.scatter(xs, [max(row['sector_leakage'], 1.0e-18) for row in sub], alpha=0.8)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')
    ax1.set_ylabel('frozen constraint')
    ax2.set_ylabel('perturbed constraint')
    ax3.set_ylabel('sector leakage')
    ax3.set_xticks(range(len(axes)), axes, rotation=25, ha='right')
    ax1.set_title('Atlas-1 constraint and leakage diagnostics')
    ax1.grid(alpha=0.25)
    ax2.grid(alpha=0.25)
    ax3.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def choose_representatives(results: list[dict[str, Any]]) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    stable = [item for item in results if item['row']['persistence_score'] == 1]
    unstable = [item for item in results if item['row']['persistence_score'] == 0]
    stable_item = None
    unstable_item = None
    if stable:
        stable_item = max(stable, key=lambda item: item['row']['coherence_score'])
    if unstable:
        unstable_item = max(unstable, key=lambda item: item['row']['baseline_coherence_score'] - item['row']['coherence_score'])
    return stable_item, unstable_item


def run_case(
    run_spec: dict[str, Any],
    baseline_row: dict[str, Any],
    defaults: dict[str, Any],
    complex_cache: dict[tuple[str, int, float], Any],
    setup_cache: dict[tuple[str, int, float, float], Any],
) -> dict[str, Any]:
    base_run = atlas_run_from_baseline(baseline_row)

    if base_run.operator_sector == 'scalar':
        key = (base_run.boundary_type, int(defaults['n_side']), float(defaults['epsilon']))
        if key not in complex_cache:
            complex_cache[key] = build_regular_complex(
                n_side=int(defaults['n_side']),
                epsilon=float(defaults['epsilon']),
                boundary_type=base_run.boundary_type,
            )
        data = complex_cache[key]
        q0, v0 = build_scalar_initial_state(
            run=base_run,
            points=data.points,
            center=np.asarray(defaults['packet_center'], dtype=float),
            kick_axis=int(defaults['kick_axis']),
            boundary_type=base_run.boundary_type,
        )
        perturbed = build_perturbed_operators(
            data=data,
            perturbation_axis=run_spec['perturbation_axis'],
            perturbation_strength=float(run_spec['perturbation_strength']),
            seed=int(run_spec['random_seed']),
            epsilon=float(defaults['epsilon']),
        )
        operator = perturbed.L0
        project = None
        divergence_op = None
        frozen_divergence_op = None
        leakage_fn = None
        positions = data.points
    else:
        key = (base_run.boundary_type, int(defaults['n_side']), float(defaults['epsilon']), float(defaults['harmonic_tol']))
        if key not in setup_cache:
            setup_cache[key] = build_transverse_setup(
                n_side=int(defaults['n_side']),
                epsilon=float(defaults['epsilon']),
                harmonic_tol=float(defaults['harmonic_tol']),
                boundary_type=base_run.boundary_type,
            )
        data, projector = setup_cache[key]
        q0, v0 = build_transverse_initial_state(
            run=base_run,
            midpoints=data.midpoints,
            edge_axes=data.edge_axes,
            center=np.asarray(defaults['packet_center'], dtype=float),
            kick_axis=int(defaults['kick_axis']),
            boundary_type=base_run.boundary_type,
            projector=projector,
        )
        perturbed = build_perturbed_operators(
            data=data,
            perturbation_axis=run_spec['perturbation_axis'],
            perturbation_strength=float(run_spec['perturbation_strength']),
            seed=int(run_spec['random_seed']),
            epsilon=float(defaults['epsilon']),
        )
        operator = perturbed.L1
        project = projector
        divergence_op = perturbed.d0.T
        frozen_divergence_op = data.d0.T
        leakage_fn = make_leakage_fn(projector)
        positions = data.midpoints

    dt = suggested_dt(operator, float(defaults['dt_safety']))
    steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
    sim = simulate_atlas_run(
        operator=operator,
        q0=q0,
        v0=v0,
        positions=positions,
        boundary_type=base_run.boundary_type,
        dt=dt,
        steps=steps,
        project=project,
        divergence_op=divergence_op,
        leakage_fn=leakage_fn,
    )

    summary = sim['summary']
    frozen_constraint_max = float(summary['max_constraint_norm'])
    perturbed_constraint_max = float(summary['max_constraint_norm'])
    if base_run.operator_sector == 'transverse':
        frozen_constraint_max = replay_constraint_max(q0, v0, dt, steps, operator, frozen_divergence_op, project=project)
        perturbed_constraint_max = replay_constraint_max(q0, v0, dt, steps, operator, divergence_op, project=project)

    labels = classify_run(
        operator_sector=base_run.operator_sector,
        summary=summary,
        max_leakage=float(summary['max_sector_leakage']),
        max_constraint=frozen_constraint_max,
        reproducible_under_refinement=False,
    )
    labels_perturbed_constraint = classify_run(
        operator_sector=base_run.operator_sector,
        summary=summary,
        max_leakage=float(summary['max_sector_leakage']),
        max_constraint=perturbed_constraint_max,
        reproducible_under_refinement=False,
    )
    perturbed_regime = labels['regime_label']
    persistence = int(perturbed_regime == run_spec['baseline_regime'])
    perturbed_regime_perturbed_constraint = labels_perturbed_constraint['regime_label']
    persistence_perturbed_constraint = int(perturbed_regime_perturbed_constraint == run_spec['baseline_regime'])

    final_q = final_values_from_sim(q0, v0, dt, steps, operator, project=project)

    row = {
        'run_id': run_spec['run_id'],
        'atlas_phase': run_spec['atlas_phase'],
        'baseline_run_id': run_spec['baseline_run_id'],
        'perturbation_axis': run_spec['perturbation_axis'],
        'perturbation_strength': float(run_spec['perturbation_strength']),
        'operator_sector': base_run.operator_sector,
        'boundary_type': base_run.boundary_type,
        'random_seed': int(run_spec['random_seed']),
        'baseline_regime': run_spec['baseline_regime'],
        'perturbed_regime': perturbed_regime,
        'perturbed_regime_perturbed_constraint': perturbed_regime_perturbed_constraint,
        'baseline_sector_label': baseline_row['sector_label'],
        'sector_label': labels['sector_label'],
        'sector_label_perturbed_constraint': labels_perturbed_constraint['sector_label'],
        'persistence_score': persistence,
        'persistence_score_perturbed_constraint': persistence_perturbed_constraint,
        'constraint_max': frozen_constraint_max,
        'constraint_max_perturbed': perturbed_constraint_max,
        'baseline_constraint_max': float(baseline_row['constraint_max']),
        'sector_leakage': float(summary['max_sector_leakage']),
        'baseline_sector_leakage': float(baseline_row['sector_leakage']),
        'norm_drift': float(summary['relative_norm_change']),
        'anisotropy_score': float(summary['max_anisotropy']),
        'coherence_score': float(summary['coherence_ratio']),
        'baseline_coherence_score': float(baseline_row['coherence_score']),
        'coherence_delta': float(summary['coherence_ratio'] - baseline_row['coherence_score']),
        'center_shift': float(summary['center_shift']),
        'width_ratio': float(summary['width_ratio']),
        'notes': run_spec['notes'],
    }
    return {
        'run': run_spec,
        'baseline': baseline_row,
        'row': row,
        'sim': sim,
        'labels': labels,
        'labels_perturbed_constraint': labels_perturbed_constraint,
        'perturbation_metadata': perturbed.metadata,
        'positions': positions,
        'final_q': final_q,
    }


def write_note(
    path: Path,
    json_rel: str,
    csv_rel: str,
    selected_specs: list[dict[str, Any]],
    regime_counts: dict[str, int],
    transition_counts: dict[str, int],
    persistence_by_regime: dict[str, float],
    axis_persistence: dict[str, float],
) -> None:
    unique_baselines = []
    for spec in selected_specs:
        if spec['baseline_run_id'] not in unique_baselines:
            unique_baselines.append(spec['baseline_run_id'])
    seed_lines = ', '.join(unique_baselines)
    path.write_text(
        '# Atlas-1 Perturbation Resilience v1\n\n'
        f'Timestamped summary: `{json_rel}`\n\n'
        f'Timestamped run table: `{csv_rel}`\n\n'
        f'Pilot baseline seeds: {seed_lines}\n\n'
        f'Perturbed regime counts: {regime_counts}\n\n'
        f'Transition counts: {transition_counts}\n\n'
        f'Persistence by baseline regime: {persistence_by_regime}\n\n'
        f'Persistence by perturbation axis: {axis_persistence}\n\n'
        'Official Atlas-1 labels are evaluated with the frozen divergence constraint because the transverse projector remains frozen from Atlas-0. '
        'The perturbed divergence constraint is retained in the run table as an auxiliary sensitivity diagnostic.\n\n'
        'Interpretation boundary: Atlas-1 maps regime persistence and transition under mild perturbation only. '
        'It does not assert continuum reconstruction, particle content, or geometry.\n',
        encoding='utf-8',
    )


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    baseline_rows = load_baseline_rows(BASELINE_CSV_PATH)
    all_runs = load_runsheet(args.runsheet)
    selected_specs = selected_runs(all_runs, args.run_ids, args.max_runs)

    complex_cache: dict[tuple[str, int, float], Any] = {}
    setup_cache: dict[tuple[str, int, float, float], Any] = {}

    results: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []

    for run_spec in selected_specs:
        baseline_row = baseline_rows[run_spec['baseline_run_id']]
        item = run_case(run_spec, baseline_row, defaults, complex_cache, setup_cache)
        results.append(item)
        csv_rows.append(item['row'])

    transition_counts = Counter(
        f"{row['baseline_regime']} -> {row['perturbed_regime']}" for row in csv_rows
    )
    regime_counts = Counter(row['perturbed_regime'] for row in csv_rows)
    persistence_by_regime = {
        label: float(np.mean([row['persistence_score'] for row in csv_rows if row['baseline_regime'] == label]))
        for label in sorted({row['baseline_regime'] for row in csv_rows})
    }
    axis_persistence = {
        axis: float(np.mean([row['persistence_score'] for row in csv_rows if row['perturbation_axis'] == axis]))
        for axis in sorted({row['perturbation_axis'] for row in csv_rows})
    }

    transition_plot = WORK_PLOT_DIR / 'stage10_atlas1_transition_heatmap.png'
    persistence_plot = WORK_PLOT_DIR / 'stage10_atlas1_persistence_by_regime.png'
    coherence_plot = WORK_PLOT_DIR / 'stage10_atlas1_coherence_change.png'
    constraint_plot = WORK_PLOT_DIR / 'stage10_atlas1_constraint_leakage.png'
    create_transition_heatmap(transition_plot, csv_rows)
    create_persistence_bar(persistence_plot, csv_rows)
    create_coherence_change_plot(coherence_plot, csv_rows)
    create_constraint_leakage_plot(constraint_plot, csv_rows)

    plot_paths = [transition_plot, persistence_plot, coherence_plot, constraint_plot]
    stable_item, unstable_item = choose_representatives(results)
    for tag, item in [('stable', stable_item), ('unstable', unstable_item)]:
        if item is None:
            continue
        run_id = item['row']['run_id']
        field_plot = WORK_PLOT_DIR / f'stage10_atlas1_{tag}_{run_id}_field_snapshot.png'
        width_plot = WORK_PLOT_DIR / f'stage10_atlas1_{tag}_{run_id}_width_trace.png'
        constraint_trace = WORK_PLOT_DIR / f'stage10_atlas1_{tag}_{run_id}_constraint_trace.png'
        create_field_snapshot(
            field_plot,
            run_id,
            item['positions'],
            item['final_q'],
            item['row']['boundary_type'],
            int(defaults['slice_axis']),
        )
        create_trace_plot(
            width_plot,
            run_id,
            item['sim']['times'],
            item['sim']['widths'],
            'packet width',
            f'Width trace ({tag})',
        )
        create_trace_plot(
            constraint_trace,
            run_id,
            item['sim']['times'],
            item['sim']['constraint_norms'],
            'constraint norm',
            f'Constraint trace ({tag})',
        )
        plot_paths.extend([field_plot, width_plot, constraint_trace])

    payload = {
        'experiment': 'stage10_atlas1_perturbation_resilience',
        'baseline_artifact': str(BASELINE_CSV_PATH.relative_to(REPO_ROOT)),
        'config': {
            'n_side': int(defaults['n_side']),
            'epsilon': float(defaults['epsilon']),
            't_final': float(defaults['t_final']),
            'dt_safety': float(defaults['dt_safety']),
            'selected_run_ids': [spec['run_id'] for spec in selected_specs],
        },
        'executed_runs': [
            {
                'run': item['run'],
                'baseline': {
                    'run_id': item['baseline']['run_id'],
                    'regime_label': item['baseline']['regime_label'],
                    'sector_label': item['baseline']['sector_label'],
                    'operator_sector': item['baseline']['operator_sector'],
                    'boundary_type': item['baseline']['boundary_type'],
                    'coherence_score': item['baseline']['coherence_score'],
                },
                'summary': item['row'],
                'labels': item['labels'],
                'labels_perturbed_constraint': item['labels_perturbed_constraint'],
                'perturbation_metadata': item['perturbation_metadata'],
                'integration': {
                    'dt': float(item['sim']['times'][1] - item['sim']['times'][0]) if len(item['sim']['times']) > 1 else 0.0,
                    'steps': max(len(item['sim']['times']) - 1, 0),
                },
            }
            for item in results
        ],
        'regime_counts': dict(regime_counts),
        'transition_counts': dict(transition_counts),
        'persistence_by_regime': persistence_by_regime,
        'persistence_by_axis': axis_persistence,
        'observation': 'Atlas-1 pilot executes cleanly on frozen Atlas-0 seeds and returns baseline-to-perturbed regime transitions under one-axis perturbations.',
        'conclusion': 'The Atlas-1 pilot is operational and can now compare regime persistence across mild perturbation channels before any threshold retuning.',
    }
    json_path, csv_path, stamped_plots, timestamp = save_atlas_payload(
        'stage10_atlas1_perturbation_resilience',
        payload,
        csv_rows,
        CSV_FIELDS,
        plot_paths,
    )

    write_note(
        NOTE_PATH,
        str(json_path.relative_to(REPO_ROOT)),
        str(csv_path.relative_to(REPO_ROOT)),
        selected_specs,
        dict(regime_counts),
        dict(transition_counts),
        persistence_by_regime,
        axis_persistence,
    )
    append_log(
        title='Stage 10B Atlas-1 Perturbation Resilience',
        config_summary=f"pilot_baseline_seeds={len({spec['baseline_run_id'] for spec in selected_specs})} perturbation_axes=edge_weight_noise,sparse_graph_defects,anisotropic_bias,kernel_width_jitter n_side={defaults['n_side']} epsilon={defaults['epsilon']}",
        result_path=json_path,
        stamped_plots=stamped_plots,
        observation='Mild one-axis perturbations shift a subset of baseline Atlas-0 seeds into neighboring morphology classes while the pipeline remains fully reproducible.',
        conclusion='Atlas-1 is operational as a local stability analysis of the morphology taxonomy and is ready for manual pilot inspection before any scaling.',
    )

    print(f'Atlas-1 pilot complete: {timestamp}')
    print(json_path)
    print(csv_path)


if __name__ == '__main__':
    main()
