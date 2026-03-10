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
    SpectrumSummary,
    append_log,
    build_dirac_system,
    degeneracy_summary,
    load_config,
    low_dirac_spectrum,
    pairing_error,
    plot_dirac_spectrum,
    timestamp_slug,
    write_json,
)


def build_note(result: dict[str, Any]) -> str:
    rows = []
    for label, record in result['systems'].items():
        rows.append(
            f"| {label} | {record['pairing_error']:.3e} | {record['lowest_abs_eigenvalue']:.3e} | {record['mean_ipr']:.3e} | {record['max_ipr']:.3e} |"
        )
    table = '\n'.join(rows)
    return f"""# DH Spectrum Structure

## Purpose

Measure the low spectral structure of the Stage 5 Dirac-type lift.

## Setup

- prototype: full chiral 4-spinor lift
- spectrum window: `{result['config']['modes']}` modes nearest zero
- systems: 2D periodic torus and 3D periodic lattice

## Spectral diagnostics

| system | pairing error | lowest | mean IPR | max IPR |
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


def summarize(summary: SpectrumSummary) -> dict[str, Any]:
    evals = np.asarray(summary.evals, dtype=float)
    return {
        'eigenvalues': evals.tolist(),
        'ipr': summary.ipr.tolist(),
        'pairing_error': pairing_error(evals),
        'lowest_abs_eigenvalue': float(np.min(np.abs(evals))),
        'mean_ipr': float(np.mean(summary.ipr)),
        'max_ipr': float(np.max(summary.ipr)),
        'degeneracy_summary': degeneracy_summary(evals),
    }


def run_dh_spectrum_structure(config: dict[str, Any] | None = None) -> tuple[dict[str, Any], Path, Path]:
    cfg = load_config(config)
    dh_cfg = cfg['DH_stage5']
    epsilon = float(cfg.get('epsilon', 0.2))
    two_d_n = int(dh_cfg.get('two_d_n', 6))
    three_d_n = int(dh_cfg.get('three_d_n', 8))
    modes = int(dh_cfg.get('spectrum_modes', 50))
    tol = float(dh_cfg.get('eig_tol', 1e-9))

    systems = {
        f'2d_periodic_n{two_d_n}': build_dirac_system(dim=2, n_side=two_d_n, epsilon=epsilon, flux_quanta=0),
        f'3d_periodic_n{three_d_n}': build_dirac_system(dim=3, n_side=three_d_n, epsilon=epsilon, flux_quanta=0),
    }
    summaries = {label: low_dirac_spectrum(system, modes=modes, tol=tol) for label, system in systems.items()}
    records = {label: summarize(summary) for label, summary in summaries.items()}

    timestamp = timestamp_slug()
    spectrum_path = PLOTS / f'{timestamp}_DH_spectrum.png'
    ipr_path = PLOTS / f'{timestamp}_DH_ipr_vs_mode.png'
    plot_dirac_spectrum(summaries, spectrum_path, ipr_path)

    observation = 'the low DH spectrum is sign-symmetric and extended in the 2D prototype, while the 3D prototype remains extended but shows only partial near-zero pairing'
    conclusion = 'the prototype shows qualitative Dirac-like spectral structure, but this does not overcome the failed square-root test against the scalar operator'
    result_path = RESULTS / f'{timestamp}_DH_spectrum_structure.json'
    result = {
        'experiment': 'DH_spectrum_structure',
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
        'plots': [str(spectrum_path.relative_to(REPO_ROOT)), str(ipr_path.relative_to(REPO_ROOT))],
        'result_path': str(result_path),
    }
    write_json(result_path, result)
    note_path = REPO_ROOT / 'experiments' / 'spinor_sector' / 'DH_Spectrum_Structure_v1.md'
    note_path.write_text(build_note(result), encoding='utf-8')
    append_log(
        title='DH Spectrum Structure',
        config_summary=f'epsilon={epsilon}, two_d_n={two_d_n}, three_d_n={three_d_n}, modes={modes}, prototype=full chiral 4-spinor lift',
        result_path=result_path,
        plot_paths=[spectrum_path, ipr_path],
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, note_path


if __name__ == '__main__':
    run_dh_spectrum_structure()
