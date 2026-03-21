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
    load_stage10_defaults,
    plt,
    save_atlas_payload,
    simulate_atlas_run,
    suggested_dt,
)
from stage10_perturbations import build_paired_perturbed_operators
from stage10_regime_labels import classify_run
from stage10_transition_labels import classify_transition_type

NOTE_PATH = ATLAS_NOTES / 'Stage_R1_Resolution_Probe_v1.md'
RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stageR1_resolution_runs.json'
BASELINE_CSV_PATH = REPO_ROOT / 'data' / '20260314_143113_stage10_atlas0_baseline_morphology.csv'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stageR1_resolution_probe_plots'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'case_id',
    'stage',
    'source_layer',
    'baseline_run_id',
    'resolution_n',
    'perturbation_mode',
    'perturbation_pair',
    'operator_sector',
    'boundary_type',
    'baseline_regime',
    'regime_label',
    'sector_label',
    'reproducible_under_refinement',
    'resolution_persistence_score',
    'transition_type_vs_base',
    'constraint_max',
    'constraint_max_aux',
    'sector_leakage',
    'norm_drift',
    'anisotropy_score',
    'coherence_score',
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
    parser = argparse.ArgumentParser(description='Run Stage R1 resolution-probing pilot.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--case-ids', nargs='*', default=[])
    parser.add_argument('--max-cases', type=int, default=0)
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


def load_runsheet(path: Path) -> tuple[list[int], list[dict[str, Any]]]:
    payload = json.loads(path.read_text(encoding='utf-8'))
    return payload['resolutions'], payload['cases']


def selected_cases(all_cases: list[dict[str, Any]], case_ids: list[str], max_cases: int) -> list[dict[str, Any]]:
    cases = all_cases
    if case_ids:
        wanted = set(case_ids)
        cases = [case for case in cases if case['case_id'] in wanted]
    if max_cases > 0:
        cases = cases[:max_cases]
    return cases


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


def create_transition_heatmap(path: Path, summaries: list[dict[str, Any]]) -> None:
    index = {label: idx for idx, label in enumerate(REGIME_ORDER)}
    matrix = np.zeros((len(REGIME_ORDER), len(REGIME_ORDER)), dtype=int)
    for item in summaries:
        matrix[index[item['base_regime']], index[item['refined_regime']]] += 1
    fig, ax = plt.subplots(figsize=(6.8, 5.5))
    im = ax.imshow(matrix, cmap='Blues')
    ax.set_xticks(range(len(REGIME_ORDER)), REGIME_ORDER, rotation=30, ha='right')
    ax.set_yticks(range(len(REGIME_ORDER)), REGIME_ORDER)
    ax.set_xlabel('refined regime (n=24)')
    ax.set_ylabel('base regime (n=16)')
    ax.set_title('Stage R1 resolution transitions')
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix[i, j]), ha='center', va='center', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_persistence_bar(path: Path, summaries: list[dict[str, Any]]) -> None:
    order = ['Ballistic coherent', 'Ballistic dispersive', 'Diffusive', 'Chaotic or irregular', 'Fragmenting']
    values = []
    for label in order:
        sub = [item['persistence_score'] for item in summaries if item['base_regime'] == label]
        values.append(float(np.mean(sub)) if sub else 0.0)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.bar(range(len(order)), values, color='tab:blue', alpha=0.82)
    ax.set_xticks(range(len(order)), order, rotation=30, ha='right')
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel('mean refinement persistence')
    ax.set_title('Stage R1 persistence by base regime')
    ax.grid(axis='y', alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_metric_delta_plot(path: Path, summaries: list[dict[str, Any]]) -> None:
    colors = {
        'stable': 'tab:blue',
        'regime_shift': 'tab:green',
        'regime_softening': 'tab:orange',
        'fragmentation_induced': 'tab:red',
        'chaos_induced': 'tab:purple',
    }
    markers = {'none': 'o', 'pair': 's'}
    fig, ax = plt.subplots(figsize=(6.3, 4.9))
    for item in summaries:
        ax.scatter(
            item['coherence_delta'],
            item['anisotropy_delta'],
            color=colors.get(item['transition_type'], 'black'),
            marker=markers.get(item['perturbation_mode'], 'o'),
            s=60,
            alpha=0.85,
        )
    ax.axvline(0.0, color='black', linestyle='--', linewidth=1.0)
    ax.axhline(0.0, color='black', linestyle='--', linewidth=1.0)
    ax.set_xlabel('coherence delta (n=24 - n=16)')
    ax.set_ylabel('anisotropy delta (n=24 - n=16)')
    ax.set_title('Stage R1 morphology drift under refinement')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_constraint_plot(path: Path, summaries: list[dict[str, Any]]) -> None:
    transverse = [item for item in summaries if item['operator_sector'] == 'transverse']
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 7.2), sharex=True)
    xs = np.arange(len(transverse))
    labels = [item['case_id'].replace('S11_', '')[:24] for item in transverse]
    ax1.plot(xs, [max(item['base_constraint'], 1.0e-18) for item in transverse], marker='o', label='frozen n=16')
    ax1.plot(xs, [max(item['refined_constraint'], 1.0e-18) for item in transverse], marker='s', label='frozen n=24')
    ax2.plot(xs, [max(item['base_constraint_aux'], 1.0e-18) for item in transverse], marker='o', label='aux n=16')
    ax2.plot(xs, [max(item['refined_constraint_aux'], 1.0e-18) for item in transverse], marker='s', label='aux n=24')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.set_ylabel('official constraint')
    ax2.set_ylabel('aux constraint')
    ax2.set_xticks(xs, labels, rotation=25, ha='right')
    ax1.set_title('Stage R1 transverse constraint scaling')
    ax1.grid(alpha=0.25)
    ax2.grid(alpha=0.25)
    ax1.legend(fontsize=8)
    ax2.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    path: Path,
    json_rel: str,
    csv_rel: str,
    selected_cases: list[dict[str, Any]],
    transition_counts: dict[str, int],
    persistence_by_regime: dict[str, float],
    persistence_by_mode: dict[str, float],
    overall_persistence: float,
) -> None:
    case_list = ', '.join(case['case_id'] for case in selected_cases)
    path.write_text(
        '# Stage R1 Resolution Probe v1\n\n'
        f'Timestamped summary: `{json_rel}`\n\n'
        f'Timestamped run table: `{csv_rel}`\n\n'
        f'Pilot cases: {case_list}\n\n'
        f'Transition counts (n=16 -> n=24): {transition_counts}\n\n'
        f'Persistence by base regime: {persistence_by_regime}\n\n'
        f'Persistence by perturbation mode: {persistence_by_mode}\n\n'
        f'Overall refinement persistence: {overall_persistence:.3f}\n\n'
        'Official transverse labels are evaluated with the frozen divergence constraint associated with the frozen projector at each resolution. '
        'Auxiliary perturbed constraints are retained only as sensitivity diagnostics for paired-stress cases.\n\n'
        'Interpretation boundary: Stage R1 probes whether atlas labels persist under refinement. '
        'It does not establish continuum limits, geometry, particles, or gauge structure.\n',
        encoding='utf-8',
    )


def run_single(
    case: dict[str, Any],
    baseline_row: dict[str, Any],
    defaults: dict[str, Any],
    n_side: int,
    complex_cache: dict[tuple[str, int, float], Any],
    setup_cache: dict[tuple[str, int, float, float], Any],
) -> dict[str, Any]:
    run = atlas_run_from_baseline(baseline_row)
    mode = case['perturbation_mode']

    if run.operator_sector == 'scalar':
        key = (run.boundary_type, n_side, float(defaults['epsilon']))
        if key not in complex_cache:
            complex_cache[key] = build_regular_complex(n_side=n_side, epsilon=float(defaults['epsilon']), boundary_type=run.boundary_type)
        data = complex_cache[key]
        q0, v0 = build_scalar_initial_state(
            run=run,
            points=data.points,
            center=np.asarray(defaults['packet_center'], dtype=float),
            kick_axis=int(defaults['kick_axis']),
            boundary_type=run.boundary_type,
        )
        if mode == 'pair':
            perturbed = build_paired_perturbed_operators(
                data=data,
                perturbations=case['component_perturbations'],
                epsilon=float(defaults['epsilon']),
                pair_name=str(case['perturbation_pair']),
            )
            operator = perturbed.L0
            aux_divergence = perturbed.d0.T
        else:
            perturbed = None
            operator = data.L0
            aux_divergence = data.d0.T
        dt = suggested_dt(operator, float(defaults['dt_safety']))
        steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
        sim = simulate_atlas_run(
            operator=operator,
            q0=q0,
            v0=v0,
            positions=data.points,
            boundary_type=run.boundary_type,
            dt=dt,
            steps=steps,
        )
        constraint_max_aux = 0.0
        max_leakage = 0.0
        max_constraint = 0.0
    else:
        key = (run.boundary_type, n_side, float(defaults['epsilon']), float(defaults['harmonic_tol']))
        if key not in setup_cache:
            setup_cache[key] = build_transverse_setup(
                n_side=n_side,
                epsilon=float(defaults['epsilon']),
                harmonic_tol=float(defaults['harmonic_tol']),
                boundary_type=run.boundary_type,
            )
        data, projector = setup_cache[key]
        q0, v0 = build_transverse_initial_state(
            run=run,
            midpoints=data.midpoints,
            edge_axes=data.edge_axes,
            center=np.asarray(defaults['packet_center'], dtype=float),
            kick_axis=int(defaults['kick_axis']),
            boundary_type=run.boundary_type,
            projector=projector,
        )
        if mode == 'pair':
            perturbed = build_paired_perturbed_operators(
                data=data,
                perturbations=case['component_perturbations'],
                epsilon=float(defaults['epsilon']),
                pair_name=str(case['perturbation_pair']),
            )
            operator = perturbed.L1
            aux_divergence = perturbed.d0.T
        else:
            perturbed = None
            operator = data.L1
            aux_divergence = data.d0.T
        dt = suggested_dt(operator, float(defaults['dt_safety']))
        steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
        sim = simulate_atlas_run(
            operator=operator,
            q0=q0,
            v0=v0,
            positions=data.midpoints,
            boundary_type=run.boundary_type,
            dt=dt,
            steps=steps,
            project=projector,
            divergence_op=data.d0.T,
            leakage_fn=make_leakage_fn(projector),
        )
        constraint_max_aux = replay_constraint_max(q0, v0, dt, steps, operator, aux_divergence, project=projector)
        max_leakage = float(sim['summary']['max_sector_leakage'])
        max_constraint = float(sim['summary']['max_constraint_norm'])

    labels = classify_run(
        operator_sector=run.operator_sector,
        summary=sim['summary'],
        max_leakage=max_leakage,
        max_constraint=max_constraint,
        reproducible_under_refinement=False,
    )

    return {
        'run': run,
        'n_side': n_side,
        'mode': mode,
        'sim': sim,
        'labels': labels,
        'constraint_max': max_constraint,
        'constraint_max_aux': float(constraint_max_aux),
        'sector_leakage': max_leakage,
        'dt': dt,
        'steps': steps,
    }


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    resolutions, all_cases = load_runsheet(args.runsheet)
    cases = selected_cases(all_cases, args.case_ids, args.max_cases)
    baseline_rows = load_baseline_rows(BASELINE_CSV_PATH)

    complex_cache: dict[tuple[str, int, float], Any] = {}
    setup_cache: dict[tuple[str, int, float, float], Any] = {}

    results_by_case: dict[str, dict[int, dict[str, Any]]] = defaultdict(dict)
    selected_case_payloads: list[dict[str, Any]] = []

    for case in cases:
        baseline_row = baseline_rows[case['baseline_run_id']]
        selected_case_payloads.append(case)
        for n_side in resolutions:
            result = run_single(case, baseline_row, defaults, int(n_side), complex_cache, setup_cache)
            results_by_case[case['case_id']][int(n_side)] = result

    run_rows: list[dict[str, Any]] = []
    case_summaries: list[dict[str, Any]] = []
    for case in cases:
        by_n = results_by_case[case['case_id']]
        base_n = int(resolutions[0])
        refined_n = int(resolutions[-1])
        base = by_n[base_n]
        refined = by_n[refined_n]

        base_regime = base['labels']['regime_label']
        refined_regime = refined['labels']['regime_label']
        persistence = int(base_regime == refined_regime)
        coherence_delta = float(refined['sim']['summary']['coherence_ratio'] - base['sim']['summary']['coherence_ratio'])
        anisotropy_delta = float(refined['sim']['summary']['max_anisotropy'] - base['sim']['summary']['max_anisotropy'])
        transition_type = classify_transition_type(
            baseline_regime=base_regime,
            perturbed_regime=refined_regime,
            baseline_width_ratio=float(base['sim']['summary']['width_ratio']),
            perturbed_width_ratio=float(refined['sim']['summary']['width_ratio']),
            coherence_delta=coherence_delta,
        )
        reproducible_under_refinement = bool(persistence)

        case_summary = {
            'case_id': case['case_id'],
            'source_layer': case['source_layer'],
            'baseline_run_id': case['baseline_run_id'],
            'operator_sector': base['run'].operator_sector,
            'boundary_type': base['run'].boundary_type,
            'perturbation_mode': case['perturbation_mode'],
            'perturbation_pair': case.get('perturbation_pair', ''),
            'base_regime': base_regime,
            'refined_regime': refined_regime,
            'base_sector_label': base['labels']['sector_label'],
            'refined_sector_label': refined['labels']['sector_label'],
            'persistence_score': persistence,
            'transition_type': transition_type,
            'coherence_delta': coherence_delta,
            'anisotropy_delta': anisotropy_delta,
            'base_constraint': float(base['constraint_max']),
            'refined_constraint': float(refined['constraint_max']),
            'base_constraint_aux': float(base['constraint_max_aux']),
            'refined_constraint_aux': float(refined['constraint_max_aux']),
            'selection_note': case['selection_note'],
        }
        case_summaries.append(case_summary)

        for n_side in resolutions:
            item = by_n[int(n_side)]
            row = {
                'run_id': f"{case['case_id']}_n{int(n_side)}",
                'case_id': case['case_id'],
                'stage': 'Stage R1',
                'source_layer': case['source_layer'],
                'baseline_run_id': case['baseline_run_id'],
                'resolution_n': int(n_side),
                'perturbation_mode': case['perturbation_mode'],
                'perturbation_pair': case.get('perturbation_pair', ''),
                'operator_sector': item['run'].operator_sector,
                'boundary_type': item['run'].boundary_type,
                'baseline_regime': base_regime,
                'regime_label': item['labels']['regime_label'],
                'sector_label': item['labels']['sector_label'],
                'reproducible_under_refinement': int(reproducible_under_refinement),
                'resolution_persistence_score': int(persistence),
                'transition_type_vs_base': 'stable' if int(n_side) == base_n else transition_type,
                'constraint_max': float(item['constraint_max']),
                'constraint_max_aux': float(item['constraint_max_aux']),
                'sector_leakage': float(item['sector_leakage']),
                'norm_drift': float(item['sim']['summary']['relative_energy_drift']),
                'anisotropy_score': float(item['sim']['summary']['max_anisotropy']),
                'coherence_score': float(item['sim']['summary']['coherence_ratio']),
                'center_shift': float(item['sim']['summary']['center_shift']),
                'width_ratio': float(item['sim']['summary']['width_ratio']),
                'notes': case['selection_note'],
            }
            run_rows.append(row)

    transition_counts = dict(Counter(f"{item['base_regime']} -> {item['refined_regime']}" for item in case_summaries))
    transition_type_counts = dict(Counter(item['transition_type'] for item in case_summaries))
    persistence_by_regime = {
        label: float(np.mean([item['persistence_score'] for item in case_summaries if item['base_regime'] == label]))
        for label in sorted({item['base_regime'] for item in case_summaries})
    }
    persistence_by_mode = {
        mode: float(np.mean([item['persistence_score'] for item in case_summaries if item['perturbation_mode'] == mode]))
        for mode in sorted({item['perturbation_mode'] for item in case_summaries})
    }
    overall_persistence = float(np.mean([item['persistence_score'] for item in case_summaries])) if case_summaries else 0.0

    plot_paths: list[Path] = []
    heatmap_path = WORK_PLOT_DIR / 'stageR1_resolution_transition_heatmap.png'
    create_transition_heatmap(heatmap_path, case_summaries)
    plot_paths.append(heatmap_path)

    bar_path = WORK_PLOT_DIR / 'stageR1_resolution_persistence_by_regime.png'
    create_persistence_bar(bar_path, case_summaries)
    plot_paths.append(bar_path)

    delta_path = WORK_PLOT_DIR / 'stageR1_resolution_metric_deltas.png'
    create_metric_delta_plot(delta_path, case_summaries)
    plot_paths.append(delta_path)

    constraint_path = WORK_PLOT_DIR / 'stageR1_resolution_constraint_scaling.png'
    create_constraint_plot(constraint_path, case_summaries)
    plot_paths.append(constraint_path)

    payload = {
        'stage': 'Stage R1',
        'description': 'Resolution-probing pilot for atlas label persistence under refinement.',
        'resolutions': resolutions,
        'cases': case_summaries,
        'transition_counts': transition_counts,
        'transition_type_counts': transition_type_counts,
        'persistence_by_regime': persistence_by_regime,
        'persistence_by_mode': persistence_by_mode,
        'overall_persistence': overall_persistence,
        'conclusion': 'Stage R1 pilot probes whether atlas labels persist from n=16 to n=24 before any broader refinement claims are attempted.',
    }

    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        experiment_slug='stageR1_resolution_probe',
        result=payload,
        csv_rows=run_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )

    json_rel = str(json_path.relative_to(REPO_ROOT))
    csv_rel = str(csv_path.relative_to(REPO_ROOT))
    write_note(
        NOTE_PATH,
        json_rel=json_rel,
        csv_rel=csv_rel,
        selected_cases=selected_case_payloads,
        transition_counts=transition_counts,
        persistence_by_regime=persistence_by_regime,
        persistence_by_mode=persistence_by_mode,
        overall_persistence=overall_persistence,
    )

    append_log(
        title=f'Stage R1 Resolution Probe ({json_path.stem})',
        config_summary=f"cases={len(cases)}, resolutions={resolutions}, csv={csv_rel}",
        result_path=json_path,
        stamped_plots=stamped_plots,
        observation=f'transition_counts={transition_counts}; overall_persistence={overall_persistence:.3f}',
        conclusion='the first Stage R1 pilot measures whether committed atlas labels survive n=16 -> n=24 refinement without changing the classifier',
    )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for plot in stamped_plots:
        print(plot)
    print(f'overall_persistence={overall_persistence:.3f}')
    print(f'transition_counts={transition_counts}')


if __name__ == '__main__':
    main()
