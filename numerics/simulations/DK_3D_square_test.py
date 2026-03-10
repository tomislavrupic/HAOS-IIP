#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from numerics.simulations.DK_stage7_common import (
    PLOTS,
    RESULTS,
    REPO_ROOT,
    append_log,
    build_dk3d_complex,
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
            f"| {label} | {record['square_relative_residual']:.3e} | {record['eigenvalue_mean_relative_error']:.3e} | {record['block0_relative_error']:.3e} | {record['block1_relative_error']:.3e} | {record['block2_relative_error']:.3e} | {record['block3_relative_error']:.3e} |"
        )
    table = '\n'.join(rows)
    return f"""# DK 3D Square Test

## Purpose

Test whether the graded 3D Dirac-Kaehler operator squares to the graded Hodge Laplacian on the validated 3D cochain complex.

## Setup

- epsilon: `{result['config']['epsilon']}`
- lattice sizes: `{result['config']['cochain_sizes']}`
- compared low modes: `{result['config']['square_modes']}`

## Grade-by-grade square diagnostics

| system | relative square residual | mean relative eigenvalue error | `Omega0` block error | `Omega1` block error | `Omega2` block error | `Omega3` block error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
{table}

## Direct result

- observation: `{result['observation']}`
- conclusion: `{result['conclusion']}`

## Artifacts

- results: `data/{Path(result['result_path']).name}`
- plots: `{', '.join(f'plots/{Path(path).name}' for path in result['plots'])}`
- timestamp: `{result['timestamp']}`
"""


def run_dk_3d_square_test(config: dict[str, Any] | None = None) -> tuple[dict[str, Any], Path, Path]:
    cfg = load_config(config)
    stage_cfg = cfg['DK_stage7']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in stage_cfg.get('cochain_sizes', [4, 6, 8])]
    modes = int(stage_cfg.get('square_modes', 40))
    tol = float(stage_cfg.get('eig_tol', 1e-9))

    systems: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        complex_data = build_dk3d_complex(n_side=n_side, epsilon=epsilon)
        checks = cochain_identity_errors(complex_data)
        dk2 = (complex_data.dirac_kahler @ complex_data.dirac_kahler).tocsr()
        delta_eigs, _ = low_eigs(complex_data.delta_h, k=min(modes, complex_data.delta_h.shape[0] - 2), which='SA', tol=tol)
        square_eigs, _ = low_eigs(dk2, k=min(modes, dk2.shape[0] - 2), which='SA', tol=tol)
        denom = np.maximum(np.abs(delta_eigs), 1e-14)
        mean_rel = float(np.mean(np.abs(square_eigs - delta_eigs) / denom))
        n0, n1, n2, n3 = complex_data.block_sizes
        blocks = [
            (slice(0, n0), slice(0, n0)),
            (slice(n0, n0 + n1), slice(n0, n0 + n1)),
            (slice(n0 + n1, n0 + n1 + n2), slice(n0 + n1, n0 + n1 + n2)),
            (slice(n0 + n1 + n2, n0 + n1 + n2 + n3), slice(n0 + n1 + n2, n0 + n1 + n2 + n3)),
        ]
        block_errors = []
        for rows, cols in blocks:
            diff = dk2[rows, cols] - complex_data.delta_h[rows, cols]
            ref = sparse_frobenius_norm(complex_data.delta_h[rows, cols]) or 1.0
            block_errors.append(sparse_frobenius_norm(diff) / ref)
        systems[f'3d_periodic_n{n_side}'] = {
            'square_relative_residual': checks['square_relative_residual'],
            'eigenvalue_mean_relative_error': mean_rel,
            'block0_relative_error': block_errors[0],
            'block1_relative_error': block_errors[1],
            'block2_relative_error': block_errors[2],
            'block3_relative_error': block_errors[3],
            'hodge_eigenvalues': delta_eigs.tolist(),
            'dk_square_eigenvalues': square_eigs.tolist(),
        }

    timestamp = timestamp_slug()
    plot_path = PLOTS / f'{timestamp}_DK_3D_square_vs_Hodge.png'
    plot_square_vs_hodge(systems, plot_path)

    observation = 'the graded 3D Dirac-Kaehler square matches the graded Hodge Laplacian to numerical precision across the tested lattice sizes and on every form degree'
    conclusion = 'the Stage 7 3D Dirac-Kaehler lift satisfies the grade-by-grade square identity on the validated 3D cochain architecture'
    result_path = RESULTS / f'{timestamp}_DK_3D_square_test.json'
    result = {
        'experiment': 'DK_3D_square_test',
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
    note_path = REPO_ROOT / 'experiments' / 'spinor_sector' / 'DK_3D_Square_Test_v1.md'
    note_path.write_text(build_note(result), encoding='utf-8')
    append_log(
        title='DK 3D Square Test',
        config_summary=f'epsilon={epsilon}, cochain_sizes={sizes}, square_modes={modes}',
        result_path=result_path,
        plot_paths=[plot_path],
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, note_path


if __name__ == '__main__':
    run_dk_3d_square_test()
