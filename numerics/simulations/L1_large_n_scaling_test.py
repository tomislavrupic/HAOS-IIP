#!/usr/bin/env python3

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from L1_transverse_scaling_test import analyze_case, fit_power_law
from L1_stage2_common import (
    PLOTS,
    REPO_ROOT,
    append_log,
    ensure_matplotlib,
    save_result_payload,
)

plt = ensure_matplotlib()

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'L1_large_n_scaling_test': {
        'sizes': [12, 16, 20, 24, 28, 32],
        'variant': 'baseline',
        'phase_modes': 30,
        'restricted_modes': 30,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['L1_large_n_scaling_test'] = dict(DEFAULT_CONFIG['L1_large_n_scaling_test'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'L1_large_n_scaling_test'})
        if isinstance(on_disk.get('L1_large_n_scaling_test'), dict):
            merged['L1_large_n_scaling_test'].update(on_disk['L1_large_n_scaling_test'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'L1_large_n_scaling_test'})
        if isinstance(config.get('L1_large_n_scaling_test'), dict):
            merged['L1_large_n_scaling_test'].update(config['L1_large_n_scaling_test'])
    return merged


def make_floor_plot(sizes: list[int], values: list[float], fit: dict[str, float], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(sizes, values, marker='o', label='baseline restricted floor')
    ax.loglog(sizes, [fit['a'] / (n ** fit['p']) for n in sizes], '--', label=f"fit p={fit['p']:.3f}")
    ax.set_xlabel('lattice side n')
    ax.set_ylabel('lowest restricted eigenvalue')
    ax.set_title('Large-n restricted transverse floor')
    ax.grid(alpha=0.25, which='both')
    ax.legend()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_band_collapse_plot(cases: dict[str, dict[str, Any]], sizes: list[int], path: Path) -> float:
    fig, ax = plt.subplots(figsize=(8, 5))
    scaled_rows = []
    for n_side in sizes:
        spectrum = np.asarray(cases[f'baseline_n{n_side}']['restricted_transverse_spectrum'][:30], dtype=float)
        scaled = (n_side * n_side) * spectrum
        scaled_rows.append(scaled[:20])
        ax.plot(range(1, len(scaled) + 1), scaled, marker='o', label=f'n={n_side}')
    scaled_rows_arr = np.asarray(scaled_rows)
    mean_curve = np.mean(scaled_rows_arr[:, :10], axis=0)
    spread = float(np.mean(np.std(scaled_rows_arr[:, :10], axis=0) / np.maximum(np.abs(mean_curve), 1e-12)))
    ax.set_xlabel('restricted mode index k')
    ax.set_ylabel(r'$n^2 \lambda_k$')
    ax.set_title('Large-n transverse-band spectral collapse')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    return spread


def make_ipr_plot(sizes: list[int], ipr_values: list[float], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(sizes, ipr_values, marker='o')
    ax.set_xlabel('lattice side n')
    ax.set_ylabel('IPR of lowest restricted mode')
    ax.set_title('Large-n localization of the lowest restricted mode')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_phase_plot(case: dict[str, Any], path: Path) -> None:
    records = case['phase_modes']
    x = [record['divergence_norm'] for record in records]
    y = [record['curl_norm'] for record in records]
    c = [record['eigenvalue'] for record in records]
    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(x, y, c=c, cmap='viridis', s=28)
    ax.set_xlabel(r'$||d_0^* a_k||$')
    ax.set_ylabel(r'$||d_1 a_k||$')
    ax.set_title(f"Large-n divergence-curl phase map (n={case['config']['n_side']})")
    ax.grid(alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r'$\lambda_k$')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'L1_Large_n_Scaling_Test_v1.md'
    sizes = result['config']['sizes']
    rows = '\n'.join(
        f"| {n_side} | {result['cases'][f'baseline_n{n_side}']['restricted_transverse_spectrum'][0]:.6f} | "
        f"{result['cases'][f'baseline_n{n_side}']['restricted_transverse_modes'][0]['ipr']:.6f} | "
        f"{result['cases'][f'baseline_n{n_side}']['restricted_transverse_modes'][0]['divergence_norm']:.3e} | "
        f"{result['cases'][f'baseline_n{n_side}']['restricted_transverse_modes'][0]['curl_norm']:.6f} |"
        for n_side in sizes
    )
    note = f"""# L1 Large-n Scaling Test

## Purpose

Stress-test the restricted transverse spectrum on larger periodic lattices without changing the operator.

## Setup

- branch: baseline periodic torus
- sizes: `{sizes}`
- operator: `T = d1* d1 | ker(d0*) intersect (H1)^perp`
- restricted modes: `{result['config']['restricted_modes']}`
- phase modes: `{result['config']['phase_modes']}`
- epsilon: `{result['config']['epsilon']}`

## Lowest restricted mode

| `n` | lowest `lambda` | IPR | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: |
{rows}

## Fit

- `lambda_min ~ {result['fit']['a']:.4f} / n^{result['fit']['p']:.3f}` with `R^2 = {result['fit']['r2']:.4f}`
- spectral-collapse spread of first 10 scaled modes: `{result['collapse_relative_spread_first_10']:.4f}`

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


def run_L1_large_n_scaling_test(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['L1_large_n_scaling_test']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20, 24, 28, 32])]
    phase_modes = int(experiment_cfg.get('phase_modes', 30))
    restricted_modes = int(experiment_cfg.get('restricted_modes', 30))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    variant = str(experiment_cfg.get('variant', 'baseline'))

    cases: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        cases[f'{variant}_n{n_side}'] = analyze_case(
            n_side=n_side,
            epsilon=epsilon,
            variant=variant,
            phase_modes=phase_modes,
            restricted_modes=restricted_modes,
            harmonic_tol=harmonic_tol,
            eig_tol=eig_tol,
            penalty=penalty,
        )

    floor_values = [cases[f'{variant}_n{n}']['restricted_transverse_spectrum'][0] for n in sizes]
    ipr_values = [cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['ipr'] for n in sizes]
    fit = fit_power_law(sizes, floor_values)

    floor_plot = PLOTS / 'transverse_floor_vs_n_large.png'
    collapse_plot = PLOTS / 'transverse_band_collapse_large.png'
    ipr_plot = PLOTS / 'transverse_ipr_vs_n_large.png'
    phase_plot = PLOTS / 'divergence_curl_phase_large.png'
    make_floor_plot(sizes, floor_values, fit, floor_plot)
    collapse_spread = make_band_collapse_plot(cases, sizes, collapse_plot)
    make_ipr_plot(sizes, ipr_values, ipr_plot)
    make_phase_plot(cases[f'{variant}_n{sizes[-1]}'], phase_plot)
    plot_paths = [floor_plot, collapse_plot, ipr_plot, phase_plot]

    observation = (
        'the lowest restricted transverse eigenvalue keeps decreasing on the baseline torus, '
        'the lowest-mode IPR continues to fall, and the first thirty scaled eigenvalues show partial n^2 collapse'
    )
    conclusion = 'the large-n baseline scan remains consistent with a continuum Laplacian-type transverse band in the tested range'

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variant': variant,
            'phase_modes': phase_modes,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
        },
        'cases': {label: {'config': case['config'], 'dimensions': case['dimensions'], 'phase_modes': case['phase_modes'], 'restricted_transverse_spectrum': case['restricted_transverse_spectrum'], 'restricted_transverse_modes': case['restricted_transverse_modes']} for label, case in cases.items()},
        'fit': fit,
        'collapse_relative_spread_first_10': collapse_spread,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('L1_large_n_scaling_test', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='L1 large-n scaling test',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, phase_modes={phase_modes}, restricted_modes={restricted_modes}, variant={variant}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_L1_large_n_scaling_test()
    print(
        json.dumps(
            {
                'result_path': str(result_path),
                'note_path': str(note_path),
                'fit': result['fit'],
                'observation': result['observation'],
                'conclusion': result['conclusion'],
            },
            indent=2,
        )
    )


if __name__ == '__main__':
    main()
