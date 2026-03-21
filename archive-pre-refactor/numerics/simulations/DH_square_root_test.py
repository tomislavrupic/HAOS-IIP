#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from numerics.simulations.DH_stage5_common import (
    PLOTS,
    RESULTS,
    REPO_ROOT,
    append_log,
    build_dirac_system,
    compare_square_root,
    load_config,
    plot_square_root_comparison,
    timestamp_slug,
    write_json,
)


def build_note(result: dict[str, Any]) -> str:
    rows = []
    for label, record in result['systems'].items():
        rows.append(
            f"| {label} | {record['mean_relative_error']:.3e} | {record['max_relative_error']:.3e} | {record['eigenvalue_correlation']:.6f} | {record['operator_relative_frobenius_error']:.3e} |"
        )
    table = '\n'.join(rows)
    return f"""# DH Square-Root Test

## Purpose

Test whether the Stage 5 Dirac-type lift behaves as a square root of the scalar operator on the same kernel substrate.

## Setup

- prototype: full chiral 4-spinor lift
- 2D baseline: `n={result['config']['two_d_n']}` periodic torus
- 3D baseline: `n={result['config']['three_d_n']}` periodic lattice
- scalar comparison: repeated low spectrum of `L0`

## Square-root diagnostics

| system | mean relative eigenvalue error | max relative eigenvalue error | eigenvalue correlation | operator relative Frobenius error |
| --- | ---: | ---: | ---: | ---: |
{table}

## Direct result

- observation: `{result['observation']}`
- conclusion: `{result['conclusion']}`

## Artifacts

- results: `data/{Path(result['result_path']).name}`
- plots: `{', '.join(f'plots/{Path(path).name}' for path in result['plots'])}`
- timestamp: `{result['timestamp']}`
"""


def run_dh_square_root_test(config: dict[str, Any] | None = None) -> tuple[dict[str, Any], Path, Path]:
    cfg = load_config(config)
    dh_cfg = cfg['DH_stage5']
    epsilon = float(cfg.get('epsilon', 0.2))
    two_d_n = int(dh_cfg.get('two_d_n', 6))
    three_d_n = int(dh_cfg.get('three_d_n', 8))
    modes = int(dh_cfg.get('square_root_modes', 48))
    tol = float(dh_cfg.get('eig_tol', 1e-9))

    systems = {
        f'2d_periodic_n{two_d_n}': build_dirac_system(dim=2, n_side=two_d_n, epsilon=epsilon, flux_quanta=0),
        f'3d_periodic_n{three_d_n}': build_dirac_system(dim=3, n_side=three_d_n, epsilon=epsilon, flux_quanta=0),
    }
    records = {label: compare_square_root(system, modes=modes, tol=tol) for label, system in systems.items()}

    timestamp = timestamp_slug()
    plot_path = PLOTS / f'{timestamp}_DH2_vs_L0_spectrum.png'
    plot_square_root_comparison(records, plot_path)

    observation = 'the chiral Dirac prototype tracks the low scalar spectrum only qualitatively: eigenvalue correlations remain high, but the relative spectral and operator errors stay large on both prototype grids'
    conclusion = 'the current Stage 5 prototype does not satisfy the square-root criterion tightly enough; under the prompt rule, this lift is rejected in its present form'
    result_path = RESULTS / f'{timestamp}_DH_square_root_test.json'
    result = {
        'experiment': 'DH_square_root_test',
        'timestamp': timestamp,
        'config': {
            'epsilon': epsilon,
            'two_d_n': two_d_n,
            'three_d_n': three_d_n,
            'modes': modes,
            'prototype': 'full chiral 4-spinor lift',
        },
        'systems': records,
        'observation': observation,
        'conclusion': conclusion,
        'plots': [str(plot_path.relative_to(REPO_ROOT))],
        'result_path': str(result_path),
    }
    write_json(result_path, result)

    note_path = Path(__file__).resolve().parents[2] / 'experiments' / 'spinor_sector' / 'DH_Square_Root_Test_v1.md'
    note_path.write_text(build_note(result), encoding='utf-8')
    append_log(
        title='DH Square-Root Test',
        config_summary=f'epsilon={epsilon}, two_d_n={two_d_n}, three_d_n={three_d_n}, modes={modes}, prototype=full chiral 4-spinor lift',
        result_path=result_path,
        plot_paths=[plot_path],
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, note_path


if __name__ == '__main__':
    run_dh_square_root_test()
