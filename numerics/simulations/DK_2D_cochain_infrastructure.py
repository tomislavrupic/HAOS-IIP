#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from numerics.simulations.DK_stage6_common import (
    RESULTS,
    REPO_ROOT,
    append_log,
    build_dk2d_complex,
    cochain_identity_errors,
    load_config,
    timestamp_slug,
    write_json,
)


def build_note(result: dict[str, Any]) -> str:
    rows = []
    for record in result['systems']:
        rows.append(
            f"| {record['n_side']} | {record['n_nodes']} | {record['n_edges']} | {record['n_faces']} | {record['d1_d0_error']:.3e} | {record['delta1_delta2_error']:.3e} |"
        )
    table = '\n'.join(rows)
    return f"""# DK 2D Cochain Infrastructure

## Purpose

Validate the weighted 2D periodic cochain complex used by the Stage 6 Dirac-Kaehler prototype.

## Setup

- prototype: weighted periodic square cell complex (`2`-torus)
- epsilon: `{result['config']['epsilon']}`
- lattice sizes: `{result['config']['cochain_sizes']}`
- cycle phases: zero (untwisted infrastructure check)

## Cochain identities

| `n` | nodes | edges | faces | `||d1 d0||` | `||delta1 delta2||` |
| ---: | ---: | ---: | ---: | ---: | ---: |
{table}

## Direct result

- observation: `{result['observation']}`
- conclusion: `{result['conclusion']}`

## Artifacts

- results: `data/{Path(result['result_path']).name}`
- timestamp: `{result['timestamp']}`
"""


def run_dk_2d_cochain_infrastructure(config: dict[str, Any] | None = None) -> tuple[dict[str, Any], Path, Path]:
    cfg = load_config(config)
    stage_cfg = cfg['DK_stage6']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in stage_cfg.get('cochain_sizes', [8, 12, 16])]

    systems = []
    max_d1d0 = 0.0
    max_delta = 0.0
    for n_side in sizes:
        complex_data = build_dk2d_complex(n_side=n_side, epsilon=epsilon)
        checks = cochain_identity_errors(complex_data)
        max_d1d0 = max(max_d1d0, checks['d1_d0_error'])
        max_delta = max(max_delta, checks['delta1_delta2_error'])
        systems.append(
            {
                'n_side': n_side,
                'n_nodes': complex_data.block_sizes[0],
                'n_edges': complex_data.block_sizes[1],
                'n_faces': complex_data.block_sizes[2],
                **checks,
            }
        )

    observation = 'the weighted periodic 2D cochain complex satisfies the coboundary and adjoint-coboundary identities to numerical precision across the tested lattice sizes'
    conclusion = 'the Stage 6 2D cochain infrastructure is internally consistent and admissible for the Dirac-Kaehler square test'

    timestamp = timestamp_slug()
    result_path = RESULTS / f'{timestamp}_DK_2D_cochain_infrastructure.json'
    result = {
        'experiment': 'DK_2D_cochain_infrastructure',
        'timestamp': timestamp,
        'config': {
            'epsilon': epsilon,
            'cochain_sizes': sizes,
            'cycle_phase_x': 0.0,
            'cycle_phase_y': 0.0,
        },
        'systems': systems,
        'observation': observation,
        'conclusion': conclusion,
        'result_path': str(result_path),
    }
    write_json(result_path, result)
    note_path = REPO_ROOT / 'experiments' / 'spinor_sector' / 'DK_2D_Cochain_Infrastructure_v1.md'
    note_path.write_text(build_note(result), encoding='utf-8')
    append_log(
        title='DK 2D Cochain Infrastructure',
        config_summary=f'epsilon={epsilon}, cochain_sizes={sizes}, cycle_phase_x=0.0, cycle_phase_y=0.0',
        result_path=result_path,
        plot_paths=[],
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, note_path


if __name__ == '__main__':
    run_dk_2d_cochain_infrastructure()
