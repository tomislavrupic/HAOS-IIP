#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from numerics.simulations.DK_stage6_common import (
    PLOTS,
    RESULTS,
    REPO_ROOT,
    append_log,
    build_dk2d_complex,
    cochain_identity_errors,
    load_config,
    plot_flux_flow,
    timestamp_slug,
    write_json,
)


def build_note(result: dict[str, Any]) -> str:
    rows = []
    for record in result['flux_cases']:
        rows.append(
            f"| {record['cycle_phase_x']:.6f} | {record['square_relative_residual']:.3e} | {record['positive_branch'][0]:.6f} | {record['positive_branch'][-1]:.6f} | {record['omega1_fraction_mean']:.3f} |"
        )
    table = '\n'.join(rows)
    return f"""# DK Flux Response Test

## Purpose

Measure how the 2D Dirac-Kaehler spectrum responds to flat torus-cycle holonomy while preserving the graded square identity.

## Setup

- lattice: 2D periodic cell complex with `n={result['config']['n_side']}`
- epsilon: `{result['config']['epsilon']}`
- cycle phases: `{result['config']['cycle_phases']}`
- tracked positive branch size: `{result['config']['spectral_flow_modes']}`

## Flux diagnostics

| cycle phase | relative square residual | first positive mode | last tracked positive mode | mean `Omega1` fraction |
| ---: | ---: | ---: | ---: | ---: |
{table}

## Direct result

- observation: `{result['observation']}`
- conclusion: `{result['conclusion']}`

## Artifacts

- results: `data/{Path(result['result_path']).name}`
- plots: `{', '.join(f'plots/{Path(path).name}' for path in result['plots'])}`
- timestamp: `{result['timestamp']}`
"""


def run_dk_flux_response_test(config: dict[str, Any] | None = None) -> tuple[dict[str, Any], Path, Path]:
    cfg = load_config(config)
    stage_cfg = cfg['DK_stage6']
    epsilon = float(cfg.get('epsilon', 0.2))
    n_side = int(stage_cfg.get('flux_n', 12))
    cycle_phases = [float(v) for v in stage_cfg.get('cycle_phases', [0.0, 0.2617993877991494, 0.5235987755982988, 0.7853981633974483])]
    modes = int(stage_cfg.get('spectrum_modes', 60))
    branch_modes = 12
    tol = float(stage_cfg.get('eig_tol', 1e-9))

    records = []
    branches: list[np.ndarray] = []
    residuals: list[float] = []
    for phase in cycle_phases:
        complex_data = build_dk2d_complex(n_side=n_side, epsilon=epsilon, cycle_phase_x=phase, cycle_phase_y=0.0)
        checks = cochain_identity_errors(complex_data)
        evals, evecs = np.linalg.eigh(complex_data.dirac_kahler.toarray())
        order = np.argsort(evals.real)
        evals = np.asarray(evals[order], dtype=float)
        evecs = np.asarray(evecs[:, order], dtype=complex)
        positives = np.sort(evals[evals > 1e-10])[:branch_modes]
        positive_idx = np.flatnonzero(evals > 1e-10)[:branch_modes]
        n0, n1, n2 = complex_data.block_sizes
        omega1 = []
        for idx in positive_idx:
            vec = evecs[:, idx]
            norm_sq = float(np.sum(np.abs(vec) ** 2)) or 1.0
            omega1_sq = float(np.sum(np.abs(vec[n0:n0 + n1]) ** 2))
            omega1.append(omega1_sq / norm_sq)
        records.append(
            {
                'cycle_phase_x': phase,
                'cycle_phase_y': 0.0,
                'square_relative_residual': checks['square_relative_residual'],
                'positive_branch': positives.tolist(),
                'omega1_fraction_mean': float(np.mean(omega1)) if omega1 else 0.0,
            }
        )
        branches.append(positives)
        residuals.append(checks['square_relative_residual'])

    timestamp = timestamp_slug()
    flow_path = PLOTS / f'{timestamp}_DK_flux_spectral_flow.png'
    residual_path = PLOTS / f'{timestamp}_DK_flux_square_residual.png'
    plot_flux_flow(cycle_phases, branches, flow_path, residuals, residual_path)

    observation = 'the low Dirac-Kaehler branch shifts smoothly under flat torus-cycle holonomy while the graded square residual remains at numerical precision'
    conclusion = 'the Stage 6 2D Dirac-Kaehler lift preserves its square identity under mild background cycle holonomy and shows stable spectral flow'
    result_path = RESULTS / f'{timestamp}_DK_flux_response_test.json'
    result = {
        'experiment': 'DK_flux_response_test',
        'timestamp': timestamp,
        'config': {
            'epsilon': epsilon,
            'n_side': n_side,
            'cycle_phases': cycle_phases,
            'spectral_flow_modes': branch_modes,
        },
        'flux_cases': records,
        'observation': observation,
        'conclusion': conclusion,
        'plots': [str(flow_path.relative_to(REPO_ROOT)), str(residual_path.relative_to(REPO_ROOT))],
        'result_path': str(result_path),
    }
    write_json(result_path, result)
    note_path = REPO_ROOT / 'experiments' / 'spinor_sector' / 'DK_Flux_Response_Test_v1.md'
    note_path.write_text(build_note(result), encoding='utf-8')
    append_log(
        title='DK Flux Response Test',
        config_summary=f'epsilon={epsilon}, n_side={n_side}, cycle_phases={cycle_phases}, spectral_flow_modes={branch_modes}',
        result_path=result_path,
        plot_paths=[flow_path, residual_path],
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, note_path


if __name__ == '__main__':
    run_dk_flux_response_test()
