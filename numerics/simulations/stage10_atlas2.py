#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from matplotlib.lines import Line2D

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
from stage10_perturbations import build_paired_perturbed_operators
from stage10_regime_labels import classify_run
from stage10_transition_labels import classify_transition_type

NOTE_PATH = ATLAS_NOTES / 'Atlas_2_Transition_Structure_v1.md'
RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage10_atlas2_runs.json'
BASELINE_CSV_PATH = REPO_ROOT / 'data' / '20260314_143113_stage10_atlas0_baseline_morphology.csv'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage10_atlas2_plots'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'atlas_phase',
    'baseline_run_id',
    'perturbation_pair',
    'operator_sector',
    'boundary_type',
    'baseline_regime',
    'perturbed_regime',
    'transition_type',
    'baseline_sector_label',
    'sector_label',
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

REGIME_ORDER = [
    'Ballistic coherent',
    'Ballistic dispersive',
    'Diffusive',
    'Localized',
    'Oscillatory trapped',
    'Chaotic or irregular',
    'Fragmenting',
    'Metastable structured',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 10C Atlas-2 transition-structure pilot.')
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
    index = {label: idx for idx, label in enumerate(REGIME_ORDER)}
    matrix = np.zeros((len(REGIME_ORDER), len(REGIME_ORDER)), dtype=int)
    for row in rows:
        matrix[index[row['baseline_regime']], index[row['perturbed_regime']]] += 1
    fig, ax = plt.subplots(figsize=(6.8, 5.5))
    im = ax.imshow(matrix, cmap='Blues')
    ax.set_xticks(range(len(REGIME_ORDER)), REGIME_ORDER, rotation=30, ha='right')
    ax.set_yticks(range(len(REGIME_ORDER)), REGIME_ORDER)
    ax.set_xlabel('perturbed regime')
    ax.set_ylabel('baseline regime')
    ax.set_title('Atlas-2 paired-stress transitions')
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix[i, j]), ha='center', va='center', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_transition_map(path: Path, rows: list[dict[str, Any]]) -> None:
    pair_markers = {
        'edge_noise_plus_anisotropic_bias': 'o',
        'defects_plus_kernel_jitter': 's',
    }
    colors = {
        'stable': 'tab:blue',
        'regime_shift': 'tab:green',
        'regime_softening': 'tab:orange',
        'fragmentation_induced': 'tab:red',
        'chaos_induced': 'tab:purple',
    }
    idx = {label: i for i, label in enumerate(REGIME_ORDER)}
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    for row in rows:
        x = idx[row['baseline_regime']]
        y = idx[row['perturbed_regime']]
        ax.scatter(
            x,
            y,
            marker=pair_markers.get(row['perturbation_pair'], 'o'),
            color=colors.get(row['transition_type'], 'black'),
            alpha=0.85,
            s=50,
        )
    ax.set_xticks(range(len(REGIME_ORDER)), REGIME_ORDER, rotation=30, ha='right')
    ax.set_yticks(range(len(REGIME_ORDER)), REGIME_ORDER)
    ax.set_xlabel('baseline regime')
    ax.set_ylabel('perturbed regime')
    ax.set_title('Atlas-2 regime transition map')
    ax.grid(alpha=0.2)

    transition_handles = [
        Line2D([0], [0], marker='o', color='w', label=label, markerfacecolor=color, markersize=8)
        for label, color in colors.items()
    ]
    pair_handles = [
        Line2D([0], [0], marker=marker, color='black', linestyle='None', label=pair, markersize=8)
        for pair, marker in pair_markers.items()
    ]
    leg1 = ax.legend(handles=transition_handles, title='transition type', fontsize=8, title_fontsize=9, loc='upper left')
    ax.add_artist(leg1)
    ax.legend(handles=pair_handles, title='pair', fontsize=8, title_fontsize=9, loc='lower right')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    path: Path,
    json_rel: str,
    csv_rel: str,
    selected_specs: list[dict[str, Any]],
    regime_counts: dict[str, int],
    transition_counts: dict[str, int],
    transition_type_counts: dict[str, int],
    persistence_by_pair: dict[str, float],
) -> None:
    unique_baselines = []
    for spec in selected_specs:
        if spec['baseline_run_id'] not in unique_baselines:
            unique_baselines.append(spec['baseline_run_id'])
    path.write_text(
        '# Atlas-2 Transition Structure v1\n\n'
        f'Timestamped summary: `{json_rel}`\n\n'
        f'Timestamped run table: `{csv_rel}`\n\n'
        f'Pilot baseline seeds: {", ".join(unique_baselines)}\n\n'
        f'Perturbed regime counts: {regime_counts}\n\n'
        f'Transition counts: {transition_counts}\n\n'
        f'Transition-type counts: {transition_type_counts}\n\n'
        f'Persistence by perturbation pair: {persistence_by_pair}\n\n'
        'Official Atlas-2 labels are evaluated with the frozen divergence constraint because the transverse projector remains frozen from Atlas-0. '
        'The paired perturbed divergence constraint is retained in the run table as an auxiliary sensitivity diagnostic.\n\n'
        'Interpretation boundary: Atlas-2 maps transition structure in regime space under paired stress only. '
        'It does not assert continuum reconstruction, particle content, or geometry.\n',
        encoding='utf-8',
    )


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
        perturbed = build_paired_perturbed_operators(
            data=data,
            perturbations=run_spec['component_perturbations'],
            epsilon=float(defaults['epsilon']),
            pair_name=run_spec['perturbation_pair'],
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
        perturbed = build_paired_perturbed_operators(
            data=data,
            perturbations=run_spec['component_perturbations'],
            epsilon=float(defaults['epsilon']),
            pair_name=run_spec['perturbation_pair'],
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

    transition_type = classify_transition_type(
        baseline_regime=run_spec['baseline_regime'],
        perturbed_regime=labels['regime_label'],
        baseline_width_ratio=float(baseline_row['width_ratio']),
        perturbed_width_ratio=float(summary['width_ratio']),
        coherence_delta=float(summary['coherence_ratio'] - baseline_row['coherence_score']),
    )

    final_q = final_values_from_sim(q0, v0, dt, steps, operator, project=project)

    row = {
        'run_id': run_spec['run_id'],
        'atlas_phase': run_spec['atlas_phase'],
        'baseline_run_id': run_spec['baseline_run_id'],
        'perturbation_pair': run_spec['perturbation_pair'],
        'operator_sector': base_run.operator_sector,
        'boundary_type': base_run.boundary_type,
        'baseline_regime': run_spec['baseline_regime'],
        'perturbed_regime': labels['regime_label'],
        'transition_type': transition_type,
        'baseline_sector_label': baseline_row['sector_label'],
        'sector_label': labels['sector_label'],
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
        'perturbation_metadata': perturbed.metadata,
        'positions': positions,
        'final_q': final_q,
    }


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
    plot_paths: list[Path] = []

    for run_spec in selected_specs:
        baseline_row = baseline_rows[run_spec['baseline_run_id']]
        item = run_case(run_spec, baseline_row, defaults, complex_cache, setup_cache)
        results.append(item)
        csv_rows.append(item['row'])

        run_id = item['row']['run_id']
        field_plot = WORK_PLOT_DIR / f'stage10_atlas2_{run_id}_field_snapshot.png'
        spectrum_plot = WORK_PLOT_DIR / f'stage10_atlas2_{run_id}_spectrum_snapshot.png'
        constraint_plot = WORK_PLOT_DIR / f'stage10_atlas2_{run_id}_constraint_trace.png'
        width_plot = WORK_PLOT_DIR / f'stage10_atlas2_{run_id}_width_trace.png'
        create_field_snapshot(field_plot, run_id, item['positions'], item['final_q'], item['row']['boundary_type'], int(defaults['slice_axis']))
        create_spectral_trace(spectrum_plot, run_id, item['sim']['times'], item['sim']['spectral_centroids'], item['sim']['spectral_spreads'])
        create_trace_plot(constraint_plot, run_id, item['sim']['times'], item['sim']['constraint_norms'], 'constraint norm', 'Constraint trace')
        create_trace_plot(width_plot, run_id, item['sim']['times'], item['sim']['widths'], 'packet width', 'Width trace')
        plot_paths.extend([field_plot, spectrum_plot, constraint_plot, width_plot])

    transition_counts = Counter(f"{row['baseline_regime']} -> {row['perturbed_regime']}" for row in csv_rows)
    regime_counts = Counter(row['perturbed_regime'] for row in csv_rows)
    transition_type_counts = Counter(row['transition_type'] for row in csv_rows)
    persistence_by_pair = {
        pair: float(np.mean([1.0 if row['transition_type'] == 'stable' else 0.0 for row in csv_rows if row['perturbation_pair'] == pair]))
        for pair in sorted({row['perturbation_pair'] for row in csv_rows})
    }

    transition_heatmap = WORK_PLOT_DIR / 'stage10_atlas2_transition_matrix.png'
    transition_map = WORK_PLOT_DIR / 'stage10_atlas2_regime_transition_map.png'
    create_transition_heatmap(transition_heatmap, csv_rows)
    create_transition_map(transition_map, csv_rows)
    plot_paths.extend([transition_heatmap, transition_map])

    payload = {
        'experiment': 'stage10_atlas2_transition_structure',
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
                    'width_ratio': item['baseline']['width_ratio'],
                },
                'summary': item['row'],
                'labels': item['labels'],
                'paired_perturbation_metadata': item['perturbation_metadata'],
                'integration': {
                    'dt': float(item['sim']['times'][1] - item['sim']['times'][0]) if len(item['sim']['times']) > 1 else 0.0,
                    'steps': max(len(item['sim']['times']) - 1, 0),
                },
            }
            for item in results
        ],
        'regime_counts': dict(regime_counts),
        'transition_counts': dict(transition_counts),
        'transition_type_counts': dict(transition_type_counts),
        'persistence_by_pair': persistence_by_pair,
        'observation': 'Atlas-2 pilot executes cleanly on frozen Atlas-1 seeds and returns regime transitions under paired perturbations.',
        'conclusion': 'The Atlas-2 pilot is operational as a small transition-structure map and can now be inspected before any scaling or Stage 11 work.',
    }
    json_path, csv_path, stamped_plots, timestamp = save_atlas_payload(
        'stage10_atlas2_transition_structure',
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
        dict(transition_type_counts),
        persistence_by_pair,
    )
    summary_plots = [
        f'plots/{timestamp}_stage10_atlas2_transition_matrix.png',
        f'plots/{timestamp}_stage10_atlas2_regime_transition_map.png',
    ]
    append_log(
        title='Stage 10C Atlas-2 Transition Structure',
        config_summary=f"pilot_baseline_seeds={len({spec['baseline_run_id'] for spec in selected_specs})} perturbation_pairs=edge_noise_plus_anisotropic_bias,defects_plus_kernel_jitter n_side={defaults['n_side']} epsilon={defaults['epsilon']}",
        result_path=json_path,
        stamped_plots=summary_plots,
        observation='Paired perturbations reveal how baseline regime labels deform under controlled stress interactions while preserving the Atlas-0/1 measurement grammar.',
        conclusion='Atlas-2 is operational as a pilot transition-topology instrument on the frozen architecture.',
    )

    print(f'Atlas-2 pilot complete: {timestamp}')
    print(json_path)
    print(csv_path)


if __name__ == '__main__':
    main()
