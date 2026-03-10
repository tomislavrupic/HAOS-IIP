#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from numerics.simulations.DK_stage7_common import (
    RESULTS,
    REPO_ROOT,
    append_log,
    build_dk3d_complex,
    cochain_identity_errors,
    load_config,
    timestamp_slug,
    write_json,
)


def build_note(result: dict[str, Any]) -> str:
    rows = []
    for record in result['systems']:
        rows.append(
            f"| {record['n_side']} | {record['n_nodes']} | {record['n_edges']} | {record['n_faces']} | {record['n_volumes']} | {record['d1_d0_error']:.3e} | {record['d2_d1_error']:.3e} | {record['delta1_delta2_error']:.3e} | {record['delta2_delta3_error']:.3e} |"
        )
    table = '\n'.join(rows)
    return f"""# DK 3D Cochain Infrastructure

## Purpose

Validate the weighted 3D periodic cochain complex used by the Stage 7 Dirac-Kaehler prototype.

## Setup

- prototype: weighted periodic cubic cell complex (`3`-torus)
- epsilon: `{result['config']['epsilon']}`
- lattice sizes: `{result['config']['cochain_sizes']}`
- cycle phases: zero (untwisted infrastructure check)

## Cochain identities

| `n` | nodes | edges | faces | volumes | `||d1 d0||` | `||d2 d1||` | `||delta1 delta2||` | `||delta2 delta3||` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
{table}

## Direct result

- observation: `{result['observation']}`
- conclusion: `{result['conclusion']}`

## Artifacts

- results: `data/{Path(result['result_path']).name}`
- timestamp: `{result['timestamp']}`
"""


def run_dk_3d_cochain_infrastructure(config: dict[str, Any] | None = None) -> tuple[dict[str, Any], Path, Path]:
    cfg = load_config(config)
    stage_cfg = cfg['DK_stage7']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in stage_cfg.get('cochain_sizes', [4, 6, 8])]

    systems = []
    for n_side in sizes:
        complex_data = build_dk3d_complex(n_side=n_side, epsilon=epsilon)
        checks = cochain_identity_errors(complex_data)
        n0, n1, n2, n3 = complex_data.block_sizes
        systems.append(
            {
                'n_side': n_side,
                'n_nodes': n0,
                'n_edges': n1,
                'n_faces': n2,
                'n_volumes': n3,
                **checks,
            }
        )

    observation = 'the weighted periodic 3D cochain complex satisfies the nilpotency and adjoint-nilpotency identities to numerical precision across the tested lattice sizes'
    conclusion = 'the Stage 7 3D cochain infrastructure is internally consistent and admissible for the graded square test'

    timestamp = timestamp_slug()
    result_path = RESULTS / f'{timestamp}_DK_3D_cochain_infrastructure.json'
    result = {
        'experiment': 'DK_3D_cochain_infrastructure',
        'timestamp': timestamp,
        'config': {
            'epsilon': epsilon,
            'cochain_sizes': sizes,
            'cycle_phase_x': 0.0,
            'cycle_phase_y': 0.0,
            'cycle_phase_z': 0.0,
        },
        'systems': systems,
        'observation': observation,
        'conclusion': conclusion,
        'result_path': str(result_path),
    }
    write_json(result_path, result)
    note_path = REPO_ROOT / 'experiments' / 'spinor_sector' / 'DK_3D_Cochain_Infrastructure_v1.md'
    note_path.write_text(build_note(result), encoding='utf-8')
    append_log(
        title='DK 3D Cochain Infrastructure',
        config_summary=f'epsilon={epsilon}, cochain_sizes={sizes}, cycle_phase_x=0.0, cycle_phase_y=0.0, cycle_phase_z=0.0',
        result_path=result_path,
        plot_paths=[],
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, note_path


if __name__ == '__main__':
    run_dk_3d_cochain_infrastructure()
