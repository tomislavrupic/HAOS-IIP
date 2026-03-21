#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import (
    ATLAS_NOTES,
    AtlasRun,
    REPO_ROOT,
    append_log,
    build_regular_complex,
    build_scalar_initial_state,
    build_transverse_initial_state,
    build_transverse_setup,
    create_field_snapshot,
    create_spectral_trace,
    create_summary_scatter,
    create_trace_plot,
    load_run_sheet,
    load_stage10_defaults,
    save_atlas_payload,
    simulate_atlas_run,
    suggested_dt,
)
from stage10_regime_labels import classify_run

NOTE_PATH = ATLAS_NOTES / 'Atlas_0_Baseline_Morphology_v1.md'
README_PATH = ATLAS_NOTES / 'README.md'
RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage10_atlas0_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage10_atlas0_plots'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'atlas_phase',
    'graph_type',
    'kernel_type',
    'operator_sector',
    'boundary_type',
    'initial_seed_type',
    'central_k',
    'bandwidth',
    'amplitude',
    'phase_pattern',
    'packet_count',
    'random_seed',
    'constraint_max',
    'sector_leakage',
    'norm_drift',
    'anisotropy_score',
    'coherence_score',
    'proto_spacetime_score',
    'regime_label',
    'sector_label',
    'notes',
    'center_shift',
    'width_ratio',
    'spectral_centroid_final',
    'spectral_spread_final',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 10 Atlas-0 baseline morphology diagnostics.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    parser.add_argument('--max-runs', type=int, default=0)
    return parser.parse_args()


def selected_runs(all_runs: list[AtlasRun], run_ids: list[str], max_runs: int) -> list[AtlasRun]:
    runs = all_runs
    if run_ids:
        wanted = set(run_ids)
        runs = [run for run in runs if run.run_id in wanted]
    if max_runs > 0:
        runs = runs[:max_runs]
    return runs


def make_leakage_fn(projector):
    def leakage(vec: np.ndarray) -> float:
        denom = max(float(np.linalg.norm(vec)), 1.0e-12)
        return float(np.linalg.norm(np.asarray(projector(vec), dtype=float) - vec) / denom)

    return leakage


def prepare_scalar_run(run: AtlasRun, defaults: dict[str, Any], complex_cache: dict[tuple[str, int, float], Any]) -> dict[str, Any]:
    key = (run.boundary_type, int(defaults['n_side']), float(defaults['epsilon']))
    if key not in complex_cache:
        complex_cache[key] = build_regular_complex(n_side=int(defaults['n_side']), epsilon=float(defaults['epsilon']), boundary_type=run.boundary_type)
    data = complex_cache[key]
    q0, v0 = build_scalar_initial_state(
        run=run,
        points=data.points,
        center=np.asarray(defaults['packet_center'], dtype=float),
        kick_axis=int(defaults['kick_axis']),
        boundary_type=run.boundary_type,
    )
    dt = suggested_dt(data.L0, float(defaults['dt_safety']))
    steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
    sim = simulate_atlas_run(
        operator=data.L0,
        q0=q0,
        v0=v0,
        positions=data.points,
        boundary_type=run.boundary_type,
        dt=dt,
        steps=steps,
    )
    return {
        'data': data,
        'sim': sim,
        'dt': dt,
        'steps': steps,
        'values_for_snapshot': np.asarray(q0 if not sim['times'] else sim['norms'], dtype=float),
    }


def prepare_transverse_run(run: AtlasRun, defaults: dict[str, Any], setup_cache: dict[tuple[str, int, float, float], Any]) -> dict[str, Any]:
    key = (run.boundary_type, int(defaults['n_side']), float(defaults['epsilon']), float(defaults['harmonic_tol']))
    if key not in setup_cache:
        setup_cache[key] = build_transverse_setup(
            n_side=int(defaults['n_side']),
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
    dt = suggested_dt(data.L1, float(defaults['dt_safety']))
    steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
    sim = simulate_atlas_run(
        operator=data.L1,
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
    return {
        'data': data,
        'sim': sim,
        'dt': dt,
        'steps': steps,
        'projector': projector,
    }


def final_values_from_sim(q0: np.ndarray, v0: np.ndarray, run_result: dict[str, Any], operator, boundary_type: str, project=None):
    q = np.asarray(q0, dtype=float).copy()
    v = np.asarray(v0, dtype=float).copy()
    if project is not None:
        q = np.asarray(project(q), dtype=float)
        v = np.asarray(project(v), dtype=float)
    dt = float(run_result['dt'])
    for _ in range(int(run_result['steps'])):
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


def run_case(run: AtlasRun, defaults: dict[str, Any], complex_cache: dict, setup_cache: dict) -> tuple[dict[str, Any], list[Path]]:
    if run.operator_sector == 'scalar':
        key = (run.boundary_type, int(defaults['n_side']), float(defaults['epsilon']))
        if key not in complex_cache:
            complex_cache[key] = build_regular_complex(n_side=int(defaults['n_side']), epsilon=float(defaults['epsilon']), boundary_type=run.boundary_type)
        data = complex_cache[key]
        q0, v0 = build_scalar_initial_state(
            run=run,
            points=data.points,
            center=np.asarray(defaults['packet_center'], dtype=float),
            kick_axis=int(defaults['kick_axis']),
            boundary_type=run.boundary_type,
        )
        dt = suggested_dt(data.L0, float(defaults['dt_safety']))
        steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
        sim = simulate_atlas_run(data.L0, q0, v0, data.points, run.boundary_type, dt, steps)
        final_q = final_values_from_sim(q0, v0, {'dt': dt, 'steps': steps}, data.L0, run.boundary_type)
        positions = data.points
        operator = data.L0
        project = None
    else:
        key = (run.boundary_type, int(defaults['n_side']), float(defaults['epsilon']), float(defaults['harmonic_tol']))
        if key not in setup_cache:
            setup_cache[key] = build_transverse_setup(
                n_side=int(defaults['n_side']),
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
        dt = suggested_dt(data.L1, float(defaults['dt_safety']))
        steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
        sim = simulate_atlas_run(
            data.L1,
            q0,
            v0,
            data.midpoints,
            run.boundary_type,
            dt,
            steps,
            project=projector,
            divergence_op=data.d0.T,
            leakage_fn=make_leakage_fn(projector),
        )
        final_q = final_values_from_sim(q0, v0, {'dt': dt, 'steps': steps}, data.L1, run.boundary_type, project=projector)
        positions = data.midpoints
        operator = data.L1
        project = projector

    summary = sim['summary']
    labels = classify_run(
        operator_sector=run.operator_sector,
        summary=summary,
        max_leakage=float(summary['max_sector_leakage']),
        max_constraint=float(summary['max_constraint_norm']),
        reproducible_under_refinement=False,
    )

    row = {
        'run_id': run.run_id,
        'atlas_phase': run.atlas_phase,
        'graph_type': run.graph_type,
        'kernel_type': run.kernel_type,
        'operator_sector': run.operator_sector,
        'boundary_type': run.boundary_type,
        'initial_seed_type': run.initial_seed_type,
        'central_k': run.central_k,
        'bandwidth': run.bandwidth,
        'amplitude': run.amplitude,
        'phase_pattern': run.phase_pattern,
        'packet_count': run.packet_count,
        'random_seed': run.random_seed,
        'constraint_max': float(summary['max_constraint_norm']),
        'sector_leakage': float(summary['max_sector_leakage']),
        'norm_drift': float(summary['relative_norm_change']),
        'anisotropy_score': float(summary['max_anisotropy']),
        'coherence_score': float(summary['coherence_ratio']),
        'proto_spacetime_score': labels['proto_spacetime_score'],
        'regime_label': labels['regime_label'],
        'sector_label': labels['sector_label'],
        'notes': run.notes,
        'center_shift': float(summary['center_shift']),
        'width_ratio': float(summary['width_ratio']),
        'spectral_centroid_final': float(summary['spectral_centroid_final']),
        'spectral_spread_final': float(summary['spectral_spread_final']),
    }

    field_plot = WORK_PLOT_DIR / f'stage10_{run.run_id}_field_snapshot.png'
    spectrum_plot = WORK_PLOT_DIR / f'stage10_{run.run_id}_spectrum_snapshot.png'
    constraint_plot = WORK_PLOT_DIR / f'stage10_{run.run_id}_constraint_trace.png'
    width_plot = WORK_PLOT_DIR / f'stage10_{run.run_id}_width_trace.png'
    create_field_snapshot(field_plot, run.run_id, positions, final_q, run.boundary_type, int(defaults['slice_axis']))
    create_spectral_trace(spectrum_plot, run.run_id, sim['times'], sim['spectral_centroids'], sim['spectral_spreads'])
    create_trace_plot(constraint_plot, run.run_id, sim['times'], sim['constraint_norms'], 'constraint norm', 'Constraint trace')
    create_trace_plot(width_plot, run.run_id, sim['times'], sim['widths'], 'packet width', 'Width trace')

    return {
        'run': run.__dict__,
        'metrics': sim,
        'summary': row,
        'labels': labels,
        'integration': {'dt': dt, 'steps': steps},
    }, [field_plot, spectrum_plot, constraint_plot, width_plot]


def write_readme() -> None:
    README_PATH.write_text(
        '# Stage 10 Pre-Geometry Atlas\n\n'
        'Atlas-0 is the baseline morphology pass for Stage 10. It reuses the frozen scalar and projected transverse operators and classifies packet behavior before interpretation.\n\n'
        'Run the full Atlas-0 sheet:\n\n'
        '`python numerics/simulations/stage10_atlas0.py`\n\n'
        'Run a small sanity subset:\n\n'
        '`python numerics/simulations/stage10_atlas0.py --max-runs 4`\n\n'
        'Outputs:\n\n'
        '- stamped JSON summary in `data/`\n'
        '- stamped CSV run table in `data/`\n'
        '- stamped plots in `plots/`\n'
        '- summary note in `experiments/pre_geometry_atlas/`\n\n'
        'Atlas-0 keeps regime labels descriptive and operator-level only. No physical interpretation is attached at this stage.\n',
        encoding='utf-8',
    )


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runs = selected_runs(load_run_sheet(args.runsheet), args.run_ids, args.max_runs)
    if not runs:
        raise SystemExit('No Atlas-0 runs selected.')

    write_readme()

    complex_cache: dict[Any, Any] = {}
    setup_cache: dict[Any, Any] = {}
    executed: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []

    for run in runs:
        result, run_plots = run_case(run, defaults, complex_cache, setup_cache)
        executed.append(result)
        csv_rows.append(result['summary'])
        plot_paths.extend(run_plots)

    summary_plot = WORK_PLOT_DIR / 'stage10_atlas0_regime_map.png'
    create_summary_scatter(summary_plot, csv_rows)
    plot_paths.append(summary_plot)

    regime_counts = Counter(row['regime_label'] for row in csv_rows)
    sector_counts = Counter(row['sector_label'] for row in csv_rows)
    observation = 'Atlas-0 executes cleanly on the frozen architecture and returns reproducible morphology labels, transport summaries, and constraint diagnostics for the selected scalar and projected transverse runs.'
    conclusion = 'The Stage 10 scaffold is in place and the Atlas-0 baseline pass can now classify packet behavior before interpretation.'
    result_payload = {
        'experiment': 'stage10_atlas0_baseline_morphology',
        'config': {
            'n_side': int(defaults['n_side']),
            'epsilon': float(defaults['epsilon']),
            't_final': float(defaults['t_final']),
            'dt_safety': float(defaults['dt_safety']),
            'selected_run_ids': [run.run_id for run in runs],
        },
        'executed_runs': executed,
        'regime_counts': dict(regime_counts),
        'sector_counts': dict(sector_counts),
        'observation': observation,
        'conclusion': conclusion,
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        experiment_slug='stage10_atlas0_baseline_morphology',
        result=result_payload,
        csv_rows=csv_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )

    append_log(
        'Stage 10 Atlas-0 Baseline Morphology',
        f"runs={len(runs)}, n_side={defaults['n_side']}, epsilon={defaults['epsilon']}, t_final={defaults['t_final']}",
        json_path,
        stamped_plots,
        observation,
        conclusion,
    )

    NOTE_PATH.write_text(
        '# Atlas-0 Baseline Morphology v1\n\n'
        f'Timestamped summary: `{json_path.relative_to(REPO_ROOT)}`\n\n'
        f'Timestamped run table: `{csv_path.relative_to(REPO_ROOT)}`\n\n'
        f'Executed runs: {", ".join(run.run_id for run in runs)}\n\n'
        f'Regime counts: {dict(regime_counts)}\n\n'
        f'Sector counts: {dict(sector_counts)}\n\n'
        'Interpretation boundary: Atlas-0 classifies packet morphology and constraint behavior only. It does not assert continuum physics or particle content.\n',
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
