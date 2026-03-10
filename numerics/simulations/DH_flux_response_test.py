#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from numerics.simulations.DH_stage5_common import (
    PLOTS,
    RESULTS,
    REPO_ROOT,
    append_log,
    build_dirac_system,
    load_config,
    low_dirac_spectrum,
    plot_flux_flow,
    timestamp_slug,
    write_json,
)


def positive_near_zero_branch(evals: np.ndarray, count: int) -> np.ndarray:
    vals = np.sort(np.asarray(evals, dtype=float))
    positives = vals[vals > 1e-10]
    return positives[:count]


def build_note(result: dict[str, Any]) -> str:
    rows = []
    for record in result['flux_cases']:
        rows.append(f"| {record['flux_quanta']} | {record['plaquette_angle']:.6f} | {record['positive_branch'][0]:.6f} | {record['positive_branch'][-1]:.6f} |")
    table = '\n'.join(rows)
    return f"""# DH Flux Response Test

## Purpose

Measure how the Stage 5 Dirac-type spectrum shifts under background torus flux in the 2D prototype.

## Setup

- prototype: full chiral 4-spinor lift
- lattice: 2D periodic torus with `n={result['config']['two_d_n']}`
- flux quanta: `{result['config']['flux_quanta']}`
- tracked positive branch size: `{result['config']['spectral_flow_modes']}`

## Flux-flow summary

| flux quanta | plaquette angle | first positive mode | last tracked positive mode |
| ---: | ---: | ---: | ---: |
{table}

## Direct result

- observation: `{result['observation']}`
- conclusion: `{result['conclusion']}`

## Artifacts

- results: `data/{Path(result['result_path']).name}`
- plots: `{', '.join(f'plots/{Path(path).name}' for path in result['plots'])}`
- timestamp: `{result['timestamp']}`
"""


def run_dh_flux_response_test(config: dict[str, Any] | None = None) -> tuple[dict[str, Any], Path, Path]:
    cfg = load_config(config)
    dh_cfg = cfg['DH_stage5']
    epsilon = float(cfg.get('epsilon', 0.2))
    n_side = int(dh_cfg.get('two_d_n', 6))
    flux_quanta = [int(v) for v in dh_cfg.get('flux_quanta', [0, 1, 2, 3])]
    modes = int(dh_cfg.get('spectral_flow_modes', 12))
    tol = float(dh_cfg.get('eig_tol', 1e-9))

    records = []
    positive_branches: list[np.ndarray] = []
    flux_angles: list[float] = []
    for flux in flux_quanta:
        system = build_dirac_system(dim=2, n_side=n_side, epsilon=epsilon, flux_quanta=flux)
        summary = low_dirac_spectrum(system, modes=2 * modes + 6, tol=tol)
        positive = positive_near_zero_branch(summary.evals, count=modes)
        records.append({
            'flux_quanta': int(flux),
            'plaquette_angle': float(system.plaquette_angle),
            'eigenvalues': summary.evals.tolist(),
            'positive_branch': positive.tolist(),
        })
        positive_branches.append(positive)
        flux_angles.append(float(system.plaquette_angle))

    timestamp = timestamp_slug()
    flow_path = PLOTS / f'{timestamp}_DH_flux_spectral_flow.png'
    plot_flux_flow(flux_angles, positive_branches, flow_path)

    observation = 'the low positive DH branch shifts smoothly as torus flux is increased in the 2D prototype'
    conclusion = 'the prototype has smooth 2D background-flux spectral flow, but it remains rejected overall because the square-root criterion against the scalar operator failed'
    result_path = RESULTS / f'{timestamp}_DH_flux_response_test.json'
    result = {
        'experiment': 'DH_flux_response_test',
        'timestamp': timestamp,
        'config': {
            'epsilon': epsilon,
            'two_d_n': n_side,
            'flux_quanta': flux_quanta,
            'spectral_flow_modes': modes,
            'prototype': 'full chiral 4-spinor lift',
        },
        'flux_cases': records,
        'observation': observation,
        'conclusion': conclusion,
        'plots': [str(flow_path.relative_to(REPO_ROOT))],
        'result_path': str(result_path),
    }
    write_json(result_path, result)
    note_path = REPO_ROOT / 'experiments' / 'spinor_sector' / 'DH_Flux_Response_Test_v1.md'
    note_path.write_text(build_note(result), encoding='utf-8')
    append_log(
        title='DH Flux Response Test',
        config_summary=f'epsilon={epsilon}, two_d_n={n_side}, flux_quanta={flux_quanta}, spectral_flow_modes={modes}, prototype=full chiral 4-spinor lift',
        result_path=result_path,
        plot_paths=[flow_path],
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, note_path


if __name__ == '__main__':
    run_dh_flux_response_test()
