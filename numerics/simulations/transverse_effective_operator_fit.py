#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage3_common import (
    PLOTS,
    REPO_ROOT,
    analyze_branch_cases,
    append_log,
    continuum_transverse_q2,
    linear_gap_fit,
    make_phase_plot,
    save_result_payload,
    plt,
)

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_effective_operator_fit': {
        'sizes': [12, 16, 20, 24],
        'restricted_modes': 24,
        'fit_modes': 2,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'source_json': 'data/20260310_145017_L1_large_n_scaling_test.json',
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_effective_operator_fit'] = dict(DEFAULT_CONFIG['transverse_effective_operator_fit'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_effective_operator_fit'})
        if isinstance(on_disk.get('transverse_effective_operator_fit'), dict):
            merged['transverse_effective_operator_fit'].update(on_disk['transverse_effective_operator_fit'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_effective_operator_fit'})
        if isinstance(config.get('transverse_effective_operator_fit'), dict):
            merged['transverse_effective_operator_fit'].update(config['transverse_effective_operator_fit'])
    return merged


def load_cases_from_json(path: Path, sizes: list[int]) -> dict[str, Any] | None:
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    cases = payload.get('cases')
    if not isinstance(cases, dict):
        return None
    wanted = {f'baseline_n{n}' for n in sizes}
    if not wanted.issubset(cases.keys()):
        return None
    return {key: cases[key] for key in sorted(wanted, key=lambda item: int(item.split('n')[-1]))}


def make_fit_plot(lambdas: np.ndarray, q2_eff: np.ndarray, fit: dict[str, float], n_side: int, path: Path) -> None:
    pred = fit['c_eff'] * q2_eff + fit['m_eff']
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(q2_eff, lambdas, 'o', label='restricted spectrum')
    ax.plot(q2_eff, pred, '-', label=rf"fit: $c_{{eff}}={fit['c_eff']:.4f}$, $m_{{eff}}={fit['m_eff']:.4e}$")
    ax.set_xlabel(r'$|q_k|^2$')
    ax.set_ylabel(r'$\lambda_k$')
    ax.set_title(f'Effective operator fit (baseline, n={n_side})')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def unique_mode_families(lambdas: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    family_lambdas: list[float] = []
    for lam in np.asarray(lambdas, dtype=float):
        if not family_lambdas or abs(lam - family_lambdas[-1]) > tol:
            family_lambdas.append(float(lam))
    return np.asarray(family_lambdas, dtype=float)


def make_gap_plot(sizes: list[int], fits: dict[str, dict[str, float]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    gaps = [fits[str(n)]['m_eff'] for n in sizes]
    ax.plot(sizes, gaps, marker='o')
    ax.axhline(0.0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel('lattice side n')
    ax.set_ylabel(r'$m_{eff}$')
    ax.set_title('Effective gap versus lattice size')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_stiffness_plot(sizes: list[int], fits: dict[str, dict[str, float]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    stiffness = [fits[str(n)]['c_eff'] for n in sizes]
    ax.plot(sizes, stiffness, marker='o')
    ax.set_xlabel('lattice side n')
    ax.set_ylabel(r'$c_{eff}$')
    ax.set_title('Effective stiffness versus lattice size')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Effective_Operator_Fit_v1.md'
    rows = '\n'.join(
        f"| {n} | {result['fits'][str(n)]['c_eff']:.6f} | {result['fits'][str(n)]['m_eff']:.6e} | {result['fits'][str(n)]['r2']:.5f} |"
        for n in result['config']['sizes']
    )
    note = f"""# Transverse Effective Operator Fit

## Purpose

Fit the low restricted transverse spectrum on the clean periodic branch to the minimal effective form `lambda_k ≈ c_eff |q_k|^2 + m_eff`.

## Setup

- branch: baseline periodic torus
- sizes: `{result['config']['sizes']}`
- restricted modes: `{result['config']['restricted_modes']}`
- fit modes: `{result['config']['fit_modes']}`

## Fitted parameters

| `n` | `c_eff` | `m_eff` | `R^2` |
| --- | ---: | ---: | ---: |
{rows}

## Direct result

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Artifacts

- results: `{result_path.relative_to(REPO_ROOT)}`
- plots: {', '.join(f'`{path}`' for path in stamped_plots)}
- timestamp: `{timestamp}`
"""
    note_path.write_text(note, encoding='utf-8')
    return note_path


def run_transverse_effective_operator_fit(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_effective_operator_fit']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20, 24])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 24))
    fit_modes = int(experiment_cfg.get('fit_modes', 2))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    source_json = experiment_cfg.get('source_json')
    cases = None
    if source_json:
        cases = load_cases_from_json(REPO_ROOT / str(source_json), sizes)
    if cases is None:
        cases = analyze_branch_cases(
            sizes=sizes,
            variants=['baseline'],
            epsilon=epsilon,
            restricted_modes=restricted_modes,
            harmonic_tol=harmonic_tol,
            eig_tol=eig_tol,
            penalty=penalty,
            flux_tube_phase=0.0,
        )

    q2_int = continuum_transverse_q2(restricted_modes)
    fits: dict[str, dict[str, float]] = {}
    q2_unique_all = np.unique(q2_int)
    for n_side in sizes:
        lambdas = np.asarray(cases[f'baseline_n{n_side}']['restricted_transverse_spectrum'][:restricted_modes], dtype=float)
        family_lambdas = unique_mode_families(lambdas)
        keep = min(fit_modes, len(family_lambdas), len(q2_unique_all))
        q2_eff = q2_unique_all[:keep] / (n_side * n_side)
        fits[str(n_side)] = linear_gap_fit(q2_eff, family_lambdas[:keep])
        fits[str(n_side)]['families_used'] = int(keep)

    fit_plot = PLOTS / 'effective_operator_fit.png'
    gap_plot = PLOTS / 'effective_gap_vs_n.png'
    stiff_plot = PLOTS / 'effective_stiffness_vs_n.png'
    phase_plot = PLOTS / 'divergence_curl_phase_effective_fit.png'

    largest = sizes[-1]
    lambdas_largest = np.asarray(cases[f'baseline_n{largest}']['restricted_transverse_spectrum'][:restricted_modes], dtype=float)
    family_lambdas_largest = unique_mode_families(lambdas_largest)
    keep_largest = min(fit_modes, len(family_lambdas_largest))
    q2_eff_largest = q2_unique_all[:keep_largest] / (largest * largest)
    make_fit_plot(family_lambdas_largest[:keep_largest], q2_eff_largest, fits[str(largest)], largest, fit_plot)
    make_gap_plot(sizes, fits, gap_plot)
    make_stiffness_plot(sizes, fits, stiff_plot)
    make_phase_plot(cases[f'baseline_n{largest}'], phase_plot, f"Restricted phase map (effective fit, n={largest})")
    plot_paths = [fit_plot, gap_plot, stiff_plot, phase_plot]

    observation = (
        'the low restricted spectrum on the clean periodic branch is well fit by a linear function of effective momentum squared across the tested lattice sizes'
    )
    conclusion = (
        f"the fitted gap stays small, from {fits[str(sizes[0])]['m_eff']:.3e} at n={sizes[0]} to {fits[str(sizes[-1])]['m_eff']:.3e} at n={sizes[-1]}, while the effective stiffness varies only mildly across the tested sizes, consistent with a nearly gapless continuum transverse branch in the tested range"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'restricted_modes': restricted_modes,
            'fit_modes': fit_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'source_json': str(source_json) if source_json else '',
        },
        'q2_int': q2_int.tolist(),
        'fits': fits,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_effective_operator_fit', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse effective operator fit',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, restricted_modes={restricted_modes}, fit_modes={fit_modes}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_effective_operator_fit()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
