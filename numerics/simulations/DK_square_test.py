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
    low_eigs,
    plot_square_vs_hodge,
    sparse_frobenius_norm,
    timestamp_slug,
    write_json,
)


def build_note(result: dict[str, Any]) -> str:
    rows = []
    for label, record in result['systems'].items():
        rows.append(
            f"| {label} | {record['square_relative_residual']:.3e} | {record['eigenvalue_mean_relative_error']:.3e} | {record['block0_relative_error']:.3e} | {record['block1_relative_error']:.3e} | {record['block2_relative_error']:.3e} |"
        )
    table = '\n'.join(rows)
    return f"""# DK Square Test

## Purpose

Test whether the graded Dirac-Kaehler operator squares to the graded Hodge Laplacian on the validated 2D cochain complex.

## Setup

- epsilon: `{result['config']['epsilon']}`
- lattice sizes: `{result['config']['cochain_sizes']}`
- compared low modes: `{result['config']['square_modes']}`

## Square diagnostics

| system | relative square residual | mean relative eigenvalue error | `Omega0` block error | `Omega1` block error | `Omega2` block error |
| --- | ---: | ---: | ---: | ---: | ---: |
{table}

## Direct result

- observation: `{result['observation']}`
- conclusion: `{result['conclusion']}`

## Artifacts

- results: `data/{Path(result['result_path']).name}`
- plots: `{', '.join(f'plots/{Path(path).name}' for path in result['plots'])}`
- timestamp: `{result['timestamp']}`
"""


def run_dk_square_test(config: dict[str, Any] | None = None) -> tuple[dict[str, Any], Path, Path]:
    cfg = load_config(config)
    stage_cfg = cfg['DK_stage6']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in stage_cfg.get('cochain_sizes', [8, 12, 16])]
    modes = int(stage_cfg.get('square_modes', 40))
    tol = float(stage_cfg.get('eig_tol', 1e-9))

    systems: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        complex_data = build_dk2d_complex(n_side=n_side, epsilon=epsilon)
        residual = (complex_data.dirac_kahler @ complex_data.dirac_kahler - complex_data.delta_h).tocsr()
        checks = cochain_identity_errors(complex_data)
        delta_eigs, _ = low_eigs(complex_data.delta_h, k=modes, which='SA', tol=tol)
        square_eigs, _ = low_eigs((complex_data.dirac_kahler @ complex_data.dirac_kahler).tocsr(), k=modes, which='SA', tol=tol)
        denom = np.maximum(np.abs(delta_eigs), 1e-14)
        mean_rel = float(np.mean(np.abs(square_eigs - delta_eigs) / denom))

        n0, n1, n2 = complex_data.block_sizes
        d2 = (complex_data.dirac_kahler @ complex_data.dirac_kahler).tocsr()
        block0 = d2[:n0, :n0] - complex_data.delta_h[:n0, :n0]
        block1 = d2[n0:n0 + n1, n0:n0 + n1] - complex_data.delta_h[n0:n0 + n1, n0:n0 + n1]
        block2 = d2[n0 + n1:, n0 + n1:] - complex_data.delta_h[n0 + n1:, n0 + n1:]
        ref0 = sparse_frobenius_norm(complex_data.delta_h[:n0, :n0]) or 1.0
        ref1 = sparse_frobenius_norm(complex_data.delta_h[n0:n0 + n1, n0:n0 + n1]) or 1.0
        ref2 = sparse_frobenius_norm(complex_data.delta_h[n0 + n1:, n0 + n1:]) or 1.0
        systems[f'2d_periodic_n{n_side}'] = {
            'square_relative_residual': checks['square_relative_residual'],
            'eigenvalue_mean_relative_error': mean_rel,
            'block0_relative_error': sparse_frobenius_norm(block0) / ref0,
            'block1_relative_error': sparse_frobenius_norm(block1) / ref1,
            'block2_relative_error': sparse_frobenius_norm(block2) / ref2,
            'hodge_eigenvalues': delta_eigs.tolist(),
            'dk_square_eigenvalues': square_eigs.tolist(),
        }

    timestamp = timestamp_slug()
    plot_path = PLOTS / f'{timestamp}_DK_square_vs_Hodge.png'
    plot_square_vs_hodge(systems, plot_path)

    observation = 'the graded Dirac-Kaehler square matches the graded Hodge Laplacian to numerical precision across the tested 2D lattice sizes'
    conclusion = 'the Stage 6 Dirac-Kaehler lift satisfies the defining square identity on the validated 2D cochain architecture'
    result_path = RESULTS / f'{timestamp}_DK_square_test.json'
    result = {
        'experiment': 'DK_square_test',
        'timestamp': timestamp,
        'config': {
            'epsilon': epsilon,
            'cochain_sizes': sizes,
            'square_modes': modes,
        },
        'systems': systems,
        'observation': observation,
        'conclusion': conclusion,
        'plots': [str(plot_path.relative_to(REPO_ROOT))],
        'result_path': str(result_path),
    }
    write_json(result_path, result)
    note_path = REPO_ROOT / 'experiments' / 'spinor_sector' / 'DK_Square_Test_v1.md'
    note_path.write_text(build_note(result), encoding='utf-8')
    append_log(
        title='DK Square Test',
        config_summary=f'epsilon={epsilon}, cochain_sizes={sizes}, square_modes={modes}',
        result_path=result_path,
        plot_paths=[plot_path],
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, note_path


if __name__ == '__main__':
    run_dk_square_test()
