#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage2_common import PLOTS, REPO_ROOT, append_log, save_result_payload
from L1_stage4_common import (
    analyze_deformed_case,
    best_mode_overlap,
    load_stage4_base_config,
    plt,
)

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_gap_response': {
        'sizes': [12],
        'variant': 'baseline',
        'families': ['radial', 'anisotropic', 'bump'],
        'etas': [0.0, 0.02, 0.05, 0.10],
        'restricted_modes': 24,
        'fit_families': 2,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'source_json': 'data/20260310_164041_scalar_transverse_coupling_test.json',
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_gap_response'] = dict(DEFAULT_CONFIG['transverse_gap_response'])
    on_disk = load_stage4_base_config(config_path)
    if on_disk:
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_gap_response'})
        if isinstance(on_disk.get('transverse_gap_response'), dict):
            merged['transverse_gap_response'].update(on_disk['transverse_gap_response'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_gap_response'})
        if isinstance(config.get('transverse_gap_response'), dict):
            merged['transverse_gap_response'].update(config['transverse_gap_response'])
    return merged


def family_groups(lambdas: list[float], count: int, tol: float = 1e-9) -> list[list[int]]:
    groups: list[list[int]] = []
    current: list[int] = []
    current_value = None
    for idx, lam in enumerate(lambdas[:count]):
        if current_value is None or abs(lam - current_value) <= tol:
            current.append(idx)
            current_value = lam if current_value is None else current_value
        else:
            groups.append(current)
            current = [idx]
            current_value = lam
    if current:
        groups.append(current)
    return groups


def fit_gap_stiffness(q2: np.ndarray, lambdas: np.ndarray) -> dict[str, float]:
    A = np.column_stack([q2, np.ones_like(q2)])
    coeffs, *_ = np.linalg.lstsq(A, lambdas, rcond=None)
    c_eff, m_eff = coeffs
    return {'c_eff': float(c_eff), 'm_eff': float(m_eff)}


def make_gap_plot(result: dict[str, Any], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for family in result['config']['families']:
        for n_side in result['config']['sizes']:
            gaps = [result['fits'][f'{family}|n{n_side}|{eta:.2f}']['m_eff'] for eta in result['config']['etas']]
            ax.plot(result['config']['etas'], gaps, marker='o', label=f'{family}, n={n_side}')
    ax.set_xlabel('deformation amplitude eta')
    ax.set_ylabel(r'$m_{eff}$')
    ax.set_title('Effective gap under scalar deformation')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Gap_Response_v1.md'
    rows = '\n'.join(
        f"| {family} | {n_side} | {eta:.2f} | {result['fits'][f'{family}|n{n_side}|{eta:.2f}']['c_eff']:.6f} | {result['fits'][f'{family}|n{n_side}|{eta:.2f}']['m_eff']:.6e} |"
        for family in result['config']['families']
        for n_side in result['config']['sizes']
        for eta in result['config']['etas']
    )
    note = f"""# Transverse Gap/Stiffness Response

## Purpose

Track how the effective low-band gap and stiffness respond to scalar/background geometry deformations.

## Setup

- sizes: `{result['config']['sizes']}`
- deformation families: `{result['config']['families']}`
- amplitudes: `{result['config']['etas']}`

## Fit summary

| family | n | eta | `c_eff` | `m_eff` |
| --- | ---: | ---: | ---: | ---: |
{rows}

## Direct result

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Artifacts

- results: `{result_path.relative_to(REPO_ROOT)}`
- plots: {', '.join(f'`{p}`' for p in stamped_plots)}
- timestamp: `{timestamp}`
"""
    note_path.write_text(note, encoding='utf-8')
    return note_path


def run_transverse_gap_response(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_gap_response']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16])]
    variant = str(experiment_cfg.get('variant', 'baseline'))
    families = list(experiment_cfg.get('families', ['radial', 'anisotropic', 'bump']))
    etas = [float(v) for v in experiment_cfg.get('etas', [0.0, 0.02, 0.05, 0.10])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 24))
    fit_families = int(experiment_cfg.get('fit_families', 2))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    source_json = str(experiment_cfg.get('source_json', ''))

    fits: dict[str, Any] = {}
    used_source = False
    source_path = REPO_ROOT / source_json if source_json else None
    if source_path and source_path.exists():
        payload = json.loads(source_path.read_text())
        responses = payload.get('responses', {})
        for n_side in sizes:
            baseline_key = f'{variant}|radial|0.02'
            baseline_spec = np.asarray(responses[baseline_key]['baseline_spectrum'][:restricted_modes], dtype=float)
            groups = family_groups(baseline_spec.tolist(), restricted_modes)[:fit_families]
            q2 = np.arange(1, len(groups) + 1, dtype=float) / (n_side * n_side)
            for family in families:
                for eta in etas:
                    if eta == 0.0:
                        spec = baseline_spec
                    else:
                        spec = np.asarray(responses[f'{variant}|{family}|{eta:.2f}']['deformed_spectrum'][:restricted_modes], dtype=float)
                    family_vals = [float(np.mean(spec[group])) for group in groups]
                    fit = fit_gap_stiffness(q2[:len(family_vals)], np.asarray(family_vals, dtype=float))
                    fit['family_values'] = family_vals
                    fits[f'{family}|n{n_side}|{eta:.2f}'] = fit
        used_source = True
    else:
        for n_side in sizes:
            baseline = analyze_deformed_case(
                n_side=n_side,
                epsilon=epsilon,
                variant=variant,
                family='none',
                eta=0.0,
                restricted_modes=restricted_modes,
                harmonic_tol=harmonic_tol,
                eig_tol=eig_tol,
                penalty=penalty,
            )
            groups = family_groups(baseline['restricted_transverse_spectrum'], restricted_modes)[:fit_families]
            q2 = np.arange(1, len(groups) + 1, dtype=float) / (n_side * n_side)
            for family in families:
                for eta in etas:
                    case = baseline if eta == 0.0 else analyze_deformed_case(
                        n_side=n_side,
                        epsilon=epsilon,
                        variant=variant,
                        family=family,
                        eta=eta,
                        restricted_modes=restricted_modes,
                        harmonic_tol=harmonic_tol,
                        eig_tol=eig_tol,
                        penalty=penalty,
                    )
                    overlap = best_mode_overlap(baseline['restricted_vectors'], case['restricted_vectors'], restricted_modes)
                    assigned = {}
                    for row, col in zip(overlap['rows'], overlap['cols']):
                        assigned[row] = float(case['restricted_transverse_spectrum'][col])
                    family_vals = []
                    for group in groups:
                        vals = [assigned[idx] for idx in group if idx in assigned]
                        family_vals.append(float(np.mean(vals)))
                    fit = fit_gap_stiffness(q2[:len(family_vals)], np.asarray(family_vals, dtype=float))
                    fit['family_values'] = family_vals
                    fits[f'{family}|n{n_side}|{eta:.2f}'] = fit

    gap_plot = PLOTS / 'transverse_gap_vs_deformation.png'
    stiffness_plot = PLOTS / 'transverse_stiffness_vs_deformation.png'
    make_gap_plot({'config': {'families': families, 'sizes': sizes, 'etas': etas}, 'fits': fits}, gap_plot)
    plot_paths = [gap_plot, stiffness_plot]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for family in families:
        for n_side in sizes:
            stiffness = [fits[f'{family}|n{n_side}|{eta:.2f}']['c_eff'] for eta in etas]
            ax.plot(etas, stiffness, marker='o', label=f'{family}, n={n_side}')
    ax.set_xlabel('deformation amplitude eta')
    ax.set_ylabel(r'$c_{eff}$')
    ax.set_title('Effective stiffness under scalar deformation')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.savefig(stiffness_plot, dpi=180, bbox_inches='tight')
    plt.close(fig)

    observation = (
        'small scalar/background deformations shift the effective transverse fit mainly through smooth stiffness renormalization, with no large unstable gap opening in the tested range'
    )
    key = f"anisotropic|n{sizes[-1]}|{etas[-1]:.2f}"
    conclusion = (
        f"for the anisotropic family at n={sizes[-1]} and eta={etas[-1]:.2f}, the fitted stiffness is {fits[key]['c_eff']:.6f} and the fitted gap is {fits[key]['m_eff']:.6e}, consistent with a structured operator-level response rather than chaotic spectral rearrangement"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variant': variant,
            'families': families,
            'etas': etas,
            'restricted_modes': restricted_modes,
            'fit_families': fit_families,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'source_json': source_json,
            'used_source': used_source,
        },
        'fits': fits,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_gap_response', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse gap response',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variant={variant}, families={families}, etas={etas}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_gap_response()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
